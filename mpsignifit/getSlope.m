function slope = getSlope(inference, cut)
% threshold = getSlope(inference, cut)
%
% Get the slope for a given cut.
%
% Inference should be a struct as generated either by BayesInference() or BootstrapInference().
%
% Cut should be the index of the desired cut.

if strcmp ( inference.call, 'bootstrap' )
    slope = inference.slopes(cut);
end;
