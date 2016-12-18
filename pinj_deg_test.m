function [KL,v,centers,edges,fitx,fitv] = pinj_deg_test(n,Pinj_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

pinj_deg = n.degree(n.p < 0);
if show_fig
    h = histogram(pinj_deg,0.5:1:max(pinj_deg)+0.5,'normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = centers;
    fitv = pdf(Pinj_dist.pd_p_neg_degree,fitx);
    hold on;
    % plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
    %     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
    plot(fitx,fitv,'-o','linewidth',2)
    legend('Histogram','Poiss.')
    xlabel('Degree')
    ylabel('Density')
    set(gca,'FontSize',20)
else
    [v,edges] = histcounts(pinj_deg,0.5:1:max(pinj_deg)+0.5,'normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = centers;
    fitv = pdf(Pinj_dist.pd_p_neg_degree,fitx);
end

%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
KL = kld(v,pdf(Pinj_dist.pd_p_neg_degree,centers),centers);
if show_fig
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f',KL);
    text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end