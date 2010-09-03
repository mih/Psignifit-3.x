function ci = getCI ( inference, cut, p )
% ci = getCI ( inference, cut )
%
% Get the confidence interval from the inference object

notin = 1-p;
probs = [0.5*notin,1-0.5*notin];

if strcmp ( inference.call, 'bootstrap' )
    bias = inference.bias_thres(cut);
    acc  = inference.acc_thres(cut);
    probs = normcdf( bias + ( norminv(probs) + bias ) ./ (1-acc*(norminv(probs) + bias )) );
end;

ci = prctile ( inference.mcthres(:,cut), 100*probs );
