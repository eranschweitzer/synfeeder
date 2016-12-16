function [KL,v,centers,edges,fitx,fitv] = voltagedrop_distribution(e,vdrop_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

vdrop_mask = e.i_est > 0 & e.inom > 0 & (e.funom==e.tunom);
vdrop = e.i_est(vdrop_mask).*sqrt(e.x(vdrop_mask).^2 + e.r(vdrop_mask).^2)./(1e3*e.tunom(vdrop_mask));
[v,edges] = histcounts(vdrop,'normalization','pdf');
centers = edges(1:end-1) + 0.5*diff(edges);
fitx = linspace(mean(edges(1:2)),edges(end),max(length(centers),100));
fitv = pdf(vdrop_dist,fitx);
KL = kld(v,pdf(vdrop_dist,centers),centers);

if show_fig
    loglog(edges(1:end-1) + 0.5*diff(edges),v,'o')
    hold on;
    plot(fitx,fitv,'linewidth',2)
    xlabel('|z|[\Omega]\times i [A]\times (u_{nom}[V])^{-1}')
    ylabel('Density')
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f',KL);
    text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
end
