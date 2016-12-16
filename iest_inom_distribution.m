function [KL,v,centers,edges,fitx,fitv] = iest_inom_distribution(e,iest_inom_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

mask = e.inom>0;
r = e.i_est(mask)./(e.inom(mask).*e.num_parallel(mask));
[v,edges] = histcounts(r,'normalization','pdf');
centers = edges(1:end-1) + 0.5*diff(edges);
fitx = linspace(edges(1),edges(end),max(length(centers),100));
fitv = pdf(iest_inom_dist,fitx);
KL = kld(v,pdf(iest_inom_dist,centers),centers);
if show_fig
    plot(centers,v,'o','MarkerFace','b')
    hold on;
    plot(fitx,fitv,'linewidth',2)
    set(gca,'Yscale','log')
    xlabel('i_{est}/i_{nom} for cables')
    ylabel('Density')
    legend('Data','Exponential distribution')
    set(gca,'FontSize',20)
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f\n Target max: 0.125',KL);
    text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
end
