function GoodnessOfFit ( inference )
% GoodnessOfFit ( inference )

ax = subplot(231);
cla(ax);
plotPMF ( inference, 'axes', ax );
ax = subplot(232);
cla(ax);
plotRd ( inference, 'p', 'axes', ax );
ax = subplot(233);
cla(ax);
plotRd ( inference, 'k', 'axes', ax );

if strcmp(inference.call, 'bayes' )
    b = inference.burnin;
    ax = subplot(234);
    cla(ax);
    hold on;
    plot ( inference.ppdeviance(b:end), inference.mcdeviance(b:end), '.' );
    yl = get ( ax, 'ylim' );
    xl = get ( ax, 'xlim' );
    ym = max(yl);
    xm = max(xl);
    plot ( [0,max(xm,ym)], [0,max(xm,ym)], 'k:' );
    xlabel ( 'simulated deviance' );
    ylabel ( 'observed deviance' );
    hold off;

    ax = subplot(235);
    cla(ax);
    hold on;
    plot ( inference.ppRpd(b:end), inference.mcRpd(b:end), '.' );
    plot ( [-1,1],[-1,1], 'k:' );
    xlabel ( 'simulated Rpd' );
    ylabel ( 'observed Rpd' );
    hold off;

    ax = subplot(236);
    cla(ax);
    hold on;
    plot ( inference.ppRkd(b:end), inference.mcRkd(b:end), '.' );
    plot ( [-1,1],[-1,1], 'k:' );
    xlabel ( 'simulated Rkd' );
    ylabel ( 'simulated Rpd' );
    hold off;
else
    diagnostics = Diagnostics ( inference.data, inference.thetahat, ...
        'sigmoid', inference.sigmoid, 'core', inference.core, ...
        'nafc', inference.nafc, 'gammaislambda', inference.gammaislambda );

    ax = subplot(234);
    cla(ax);
    hold on;
    hist(inference.mcdeviance);
    yl = get ( ax, 'ylim' );
    p = prctile ( inference.mcdeviance(:), 95 );
    plot ( [p,p], yl, 'r:' );
    plot ( [diagnostics.deviance,diagnostics.deviance], yl, 'r', 'linewidth', 2 );
    xlabel ( 'deviance' );
    ylabel ( 'number of occuances' );
    hold off;

    ax = subplot(235);
    cla(ax);
    hold on;
    hist ( inference.mcRpd );
    yl = get ( ax, 'ylim' );
    p = prctile ( inference.mcRpd(:), [2.5,97.5] );
    plot ( [p(1),p(1)], yl, 'r:', [p(2),p(2)], yl, 'r:' );
    plot ( [diagnostics.rpd,diagnostics.rpd], yl, 'r', 'linewidth', 2 );
    xlabel ( 'Rpd' );
    ylabel ( 'number of occuances' );
    hold off;

    ax = subplot(236);
    cla(ax);
    hold on;
    hist ( inference.mcRkd );
    yl = get ( ax, 'ylim' );
    p = prctile ( inference.mcRkd(:), [2.5,97.5] );
    plot ( [p(1),p(1)], yl, 'r:', [p(2),p(2)], yl, 'r:' );
    plot ( [diagnostics.rkd,diagnostics.rkd], yl, 'r', 'linewidth', 2 );
    xlabel ( 'Rkd' );
    ylabel ( 'number of occuances' );
    hold off;
end
