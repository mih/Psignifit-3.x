function results = MapEstimate ( data, priors, varargin )
% results = MapEstimate ( data, priors, ... )
%
% See BootstrapInference for additional parameters
%
% This function is part of psignifit3 (c) 2010 by Ingo Fründ

psignifitpath = '../cli';

% Check data format
if size ( data )(2) != 3
    error ( 'data should have three columns' );
end

% default values
nafc = 2;
sigmoid = 'logistic';
core    = 'mw0.1';
gil = '';
gammaislambda = false;
verbosity = '';
cuts = [0.25,0.5,0.75];
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
    case 'verbose'
        verbosity = '-v';
        verbose = true;
    case 'cuts'
        [cuts,varargin] = popoption(varargin);
    otherwise
        printf ( 'unknown option: %s !\n' , char(opt) );
    end
end

% Store the data
save ( '-ascii', '__data.txt', 'data' );

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
cmd = sprintf ( '%s/psignifit-mapestimate __data.txt --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -s %s -c %s -cuts "%s"', ...
    psignifitpath,
    getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    sigmoid, core, scuts);

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

results.call = 'mapestimate';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.burnin = 1;
results.nsamples = 0;

% clean up
% delete ( '__data.txt' );
% delete ( '__mapestimate.m' );
