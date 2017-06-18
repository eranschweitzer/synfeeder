function [KL,v,centers,edges,fitx,fitv] = pinj_frac_test(n,Pinj_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
pinj_frac = zeros(length(unique(n.fid)),1);
for fid = unique(n.fid).'
    fid_mask = n.fid == fid;
    pinj_frac(fid) = sum((n.p < 0) & fid_mask)/sum(fid_mask);
end
if show_fig
    h = histogram(pinj_frac,'normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end),max(length(v),100));
    fitv = pdf(Pinj_dist.pd_p_neg_frac,fitx);
    hold on;
    % plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
    %     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
    plot(fitx,fitv,'linewidth',2)
    legend('Histogram','Beta fit')
    xlabel('Fraction of Injection Nodes')
    ylabel('Density')
    set(gca,'FontSize',20)
else
    [v,edges] = histcounts(pinj_frac,'normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end),max(length(v),100));
    fitv = pdf(Pinj_dist.pd_p_neg_frac,fitx);
end

%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
KL = kld(v,pdf(Pinj_dist.pd_p_neg_frac,centers),centers);
if show_fig
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f',KL);
    text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end