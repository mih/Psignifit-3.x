function results = BootstrapInference ( data, priors, varargin )
% results = BootstrapInference ( data, priors, ... )
%
% Possible additional parameters
% ------------------------------
%
% 'nafc', integer       Default: 2
%     the keyword 'nafc' followed by an integer can be used to specify a yes-no task.
% 'samples', integer    Default: 2000
%     the keyword 'samples' followed by an integer can be used to specify the number of samples
% 'sigmoid', string     Default: 'logistic'
%     the keyword 'sigmoid' can be used to specify the form of the sigmoidal function. Valid choices are
%     'logistic', 'gauss', 'gumbel_l', 'gumbel_r', 'exponential', and 'cauchy'
% 'core', string        Default: 'mw0.1'
%     the keyword 'core' can be used to specify the core object of the psychometric function. Valid choices
%     are 'ab', 'mw', 'linear', 'log', 'poly', 'weibull'
% 'cuts', vector        Default: [0.25,0.5,0.75]
%     use this to specify the cuts, i.e. the positions at which thresholds should be determined
% 'gammaislambda'
%     switch on the constraint that gamma should always equal lambda
% 'verbose'
%     switch on higher verbosity
% 'nonparametric'
%     perform nonparametric bootstrap instead of the default parametric bootstrap
%
% Examples of usage:
% ------------------
%
% results = BootstrapInference ( data, priors, 'nafc', 1, 'gammaislambda' )
%    fits a yes-no task, constraining gamma to equal lambda
%
% results = BootstrapInference ( data, priors, 'nonparametric' )
%    performs nonparametric bootstrap
%
% This function is part of psignifit3 (c) 2010 by Ingo Fründ

psignifitpath = '../cli';

% Check data format
if size ( data, 2 ) ~= 3
    error ( 'data should have three columns' );
end

% default values
nafc = 2;
sigmoid = 'logistic';
core    = 'mw0.1';
gil = '';
gammaislambda = false;
npr = '';
verbosity = '';
cuts = [0.25,0.5,0.75];
samples = 2000;
verbose = false;

% Check input
while size(varargin,2) > 0
    [opt,varargin] = popoption ( varargin );
    switch opt
    case 'nafc'
        [nafc,varargin] = popoption(varargin);
    case 'sigmoid'
        [sigmoid,varargin] = popoption(varargin);
    case 'core'
        [core,varargin] = popoption(varargin);
    case 'gammaislambda'
        gil = '-gammaislambda';
        gammaislambda = true;
    case 'nonparametric'
        npr = '-nonparametric';
    case 'verbose'
        verbosity = '-v';
        verbose = true;
    case 'cuts'
        [cuts,varargin] = popoption(varargin);
    case 'samples'
        [samples,varargin] = popoption(varargin);
    otherwise
        printf ( 'unknown option: %s !\n' , char(opt) );
    end
end

% Get the point estimate
if gammaislambda
    mapest = MapEstimate ( data, priors, 'nafc', nafc, 'sigmoid', sigmoid, 'core', core, 'gammaislambda' );
else
    mapest = MapEstimate ( data, priors, 'nafc', nafc, 'sigmoid', sigmoid, 'core', core );
end

% Store the data
save ( '-ascii', '__data.txt', 'data' );

% Fiddle around with the fourth prior. Do we need it?
if nafc > 1
    prior4 = '';
elseif gammaislambda
    prior4 = '';
else
    prior4 = sprintf ( '-prior4 "%s"', getfield ( priors, 'gamma' ) );
end

% Determine cuts
scuts = sprintf ( '"%s', num2str ( cuts, '%f,') );
scuts(length(scuts)) = '"';

% Write the command
cmd = sprintf ( '%s/psignifit-bootstrap %s __data.txt --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -nsamples %d -s %s -c %s %s %s -cuts %s', ...
    psignifitpath, verbosity, ...
    getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    samples, sigmoid, core, gil, npr, scuts );

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

% Store paradigm
results.call = 'bootstrap';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.thetahat = mapest.thetahat;
results.mapest = mapest;
results.burnin = 1;
results.nsamples = samples;

% Clean up
delete ( '__data.txt' );
