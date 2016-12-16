function [KL,v,centers,edges,fitx,fitv] = hop_distance_distribution(n,hop_dist,varargin)

% plot(x_d_hop,d_hist.hop.values,'o','MarkerFaceColor','b');
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
if show_fig
    h = histogram(n.d_hop,-0.5:1:max(n.d_hop)+0.5,'Normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    hold on;
    % plot(x_d_hop,nbinpdf(x_d_hop,phat(1),phat(2)),'linewidth',2);
    plot(centers,pdf(hop_dist,centers),'linewidth',2);
    xlabel('Distance To Source [hop]')
    ylabel('Density')
    legend('Histogram','Negative Binomial Fit')
    set(gca,'FontSize',20)
    set(gca,'Xlim',[0 max(n.d_hop)])
else
    [v,edges]=histcounts(n.d_hop,-0.5:1:max(n.d_hop)+0.5,'Normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
end

%% KL Divergence of hop distance
KL = kld(v,pdf(hop_dist,centers),centers);
ax = gca;
if show_fig
    str = sprintf('KL-Divergence: %0.4f\n <90%%: 0.39',KL);
    text(ax.XLim(1),ax.YLim(2),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','left')
end


fitx = centers.';
fitv = pdf(hop_dist,centers).';