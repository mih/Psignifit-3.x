function results = BayesInference ( data, priors, varargin )
% results = BayesInference ( data, priors, sample=2000, cuts=[0.25, 0.5, 0.75], core='mw0.1', sigmoid='logistic', nafc=2
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
stepwidths = [0.1,0.1,0.01];
pilot = false;
generic = '';

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
    case 'stepwidths'
        [stepwidths,varargin] = popoption(varargin);
    case 'pilot'
        pilot = true;
        [pilotsample,varargin] = popoption(varargin);
        f = fopen ( '__pilot.txt', 'w' );
        fprintf ( f, '\n# mcestimates\n' );
        for k = 1:size(pilotsample,1)
            fprintf ( f, '%s\n', num2str(pilotsample(k,:)) );
        end
        fclose(f);
    case 'generic'
        generic = '-generic';
    otherwise
        printf ( 'unknown option: %s !\n' , char(opt) );
    end
end

if pilot
    stepwidths_or_pilot = '__pilot.txt';
else
    stepwidths_or_pilot = sprintf('"%s',num2str(stepwidths(:)','%f,'));
    stepwidths_or_pilot(end) = '"';
end

stepwidths_or_pilot

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
cmd = sprintf ( '%s/psignifit-mcmc %s __data.txt --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -nsamples %d -s %s -c %s %s %s -cuts %s -proposal %s %s', ...
    psignifitpath, verbosity, ...
    getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    samples, sigmoid, core, gil, npr, scuts, stepwidths_or_pilot, generic );

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

% Store paradigm
results.call = 'bayes';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.thetahat = mean ( results.mcestimates(samples/2:end,:) );
results.burnin = samples/2;
results.nsamples = samples;

% Clean up
delete ( '__data.txt' );
