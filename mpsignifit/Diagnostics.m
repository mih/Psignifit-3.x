function results = Diagnostics ( data, parameters, varargin )
% results = Diagnostics ( data, parameters, ... )
%
% Determine some parameters of fitted psychometric function parameters
% See BootstrapInference for more information about additional parameters
%
% This file is 

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

sparams = sprintf ( '"%s', num2str ( parameters, '%f,') );
sparams(end) = '"';

% Determine cuts
scuts = sprintf ( '"%s', num2str ( cuts, '%f,') );
scuts(end) = '"';

% Write the command
cmd = sprintf ( '%s/psignifit-diagnostics __data.txt --matlab -c %s -s %s -params %s -cuts %s -nafc %d %s', ...
    psignifitpath, core, sigmoid,sparams,scuts,nafc,verbosity );

if verbose
    cmd
end

[status,output] = system ( cmd, 1 );
eval ( output );

results.call = 'diagnostics';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;

% Clean up
% delete ( '__data.txt' );
% delete ( sprintf('%s.m',name) );
