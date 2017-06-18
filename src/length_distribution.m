function [KL,v,centers,edges,fitx,fitv] = length_distribution(e,length_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
[v,edges] = histcounts(e.length(e.length>0),50,'normalization','pdf');
centers = edges(1:end-1) + 0.5*diff(edges);
fitx = linspace(mean(edges(1:2)),edges(end),max(length(centers),100));
fitv = length_dist.cauchyMod(fitx,length_dist.cauchyHat(1),length_dist.cauchyHat(2));
KL = kld(v,length_dist.cauchyMod(centers,length_dist.cauchyHat(1),length_dist.cauchyHat(2)),centers);

if show_fig
    loglog(edges(1:end-1) + 0.5*diff(edges),v,'o')
    hold on; 
    plot(linspace(mean(edges(1:2)),edges(end)),...
        length_dist.cauchyMod(linspace(mean(edges(1:2)),edges(end)),...
        length_dist.cauchyHat(1),length_dist.cauchyHat(2)))
    xlabel('length [km]')
    ylabel('Density')
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f\n Target max: 0.3',KL);
    text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
end
