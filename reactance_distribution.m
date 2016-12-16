function [KL,v,centers,edges,fitx,fitv] = reactance_distribution(e,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

cauchyMod = @(data,x0,gamma) (atan(x0/gamma) + pi/2)^(-1).*(gamma./((data-x0).^2 + gamma.^2));
% Fit modified cauchy distribution to data
% cauchyHat = mle(e.x(e.x>0),'pdf',cauchyMod,'start',[0 1],'lower',[-inf eps],'upper',[inf inf]);
cauchyHat = [0.0446,0.0351];


[v,edges] = histcounts(e.x(e.x>0),50,'normalization','pdf');
centers = edges(1:end-1) + 0.5*diff(edges);
fitx = linspace(mean(edges(1:2)),edges(end),max(length(centers),100));
fitv = cauchyMod(fitx,cauchyHat(1),cauchyHat(2));
KL = kld(v,cauchyMod(centers,cauchyHat(1),cauchyHat(2)),centers);

if show_fig
    loglog(edges(1:end-1) + 0.5*diff(edges),v,'o')

    hold on; 
    plot(linspace(mean(edges(1:2)),edges(end)),...
        cauchyMod(linspace(mean(edges(1:2)),edges(end)),cauchyHat(1),cauchyHat(2)))
    xlabel('Reactance [\Omega]')
    ylabel('Density')
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f\n Target max: ',KL);
    text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
end