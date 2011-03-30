function threshold = getThres(inference, cut)
% threshold = getThreshold(inference, cut)
%
% Get the threshold for a given cut.
%
% Inference should be a struct as generated either by BayesInference() or BootstrapInference().
%
% Cut should be the index of the desired cut.

if strcmp ( inference.call, 'bootstrap' )
    threshold = inference.thresholds(cut);
end;
