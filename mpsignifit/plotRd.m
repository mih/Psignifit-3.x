function axhandle = plotRd ( inference, regressor, varargin )
% axhandle = plotRd ( inference, regressor, varargin )
%
%

% Check data format
if size ( inference.data )(2) != 3
    error ( 'data should have three columns' );
end

if !isstruct ( inference )
    error ( 'inference should be a struct' );
end

axhandle = gca;
color = 'b';
if strcmp(regressor, 'p')
    xname = 'model prediction';
elseif strcmp(regressor, 'k')
    xame = 'block index';
else
    error ( sprintf ( 'Unknown regressor %s', regressor ) );
end
xname = 'block index';
yname = 'deviance residual';

while size(varargin,2) > 0
    [opt,varargin] = popoption ( varargin );
    switch opt
    case 'verbose'
        verbose = true;
    case 'axes'
        [axhandle,varargin] = popoption ( varargin );
    case 'color'
        [color,varargin] = popoption ( varargin );
    case 'xlabel'
        [xname,varargin] = popoption ( varargin );
    case 'ylabel'
        [yname,varargin] = popoption ( varargin );
    otherwise
        printf ( 'unknown option: %s !\n' , char(opt) );
    end
end


diagnostics = Diagnostics ( inference.data, inference.thetahat, ...
    'sigmoid', inference.sigmoid, 'core', inference.core, ...
    'nafc', inference.nafc, 'gammaislambda', inference.gammaislambda );

if strcmp(regressor, 'p')
    x = diagnostics.prediction(:,2);
    R = diagnostics.rpd;
elseif strcmp(regressor, 'k' )
    x = 1:size(diagnostics.prediction,1);
    R = diagnostics.rkd;
else
    error ( sprintf ( 'Unknown regressor %s', regressor ) );
end

b = cov ( x, diagnostics.devianceresiduals ) ./ var ( x );
a = mean(diagnostics.devianceresiduals) - b*mean(x);

r = cov ( x, diagnostics.devianceresiduals ) ./ sqrt( var(x).*var(diagnostics.devianceresiduals) );

cla(axhandle);

hold on;
plot ( axhandle, x, diagnostics.devianceresiduals, 'o', 'color', color );
x = linspace ( min(x), max(x) );
plot ( axhandle, x, a + b*x, ':', 'color', color );
hold off;

d = 0.05*(max(x)-min(x));
set ( axhandle, 'xlim', [ min(x)-d, max(x)+d ] );

xlabel ( xname );
ylabel ( yname );
text ( min(x)+18*d, min(diagnostics.devianceresiduals) + 0.9*(max(diagnostics.devianceresiduals)-min(diagnostics.devianceresiduals)), ...
    sprintf ( 'R%sd=%.2f', regressor, R) );
