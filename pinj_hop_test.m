function [KL,v,centers,edges,fitx,fitv] = pinj_hop_test(n,Pinj_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

pinj_hopfrac = [];
for fid = unique(n.fid).'
    fid_mask = n.fid == fid;
    Max_hop = max(n.d_hop(fid_mask));
    pinj_mask = (n.p < 0) & fid_mask;
    pinj_hopfrac = [pinj_hopfrac; (n.d_hop(pinj_mask)./Max_hop)]; %#ok<AGROW>
end
if show_fig
    h = histogram(pinj_hopfrac,'normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end),max(length(v),100));
    fitv = Pinj_dist.norm_bi(fitx,Pinj_dist.phat(1),Pinj_dist.phat(2),...
            Pinj_dist.phat(3),Pinj_dist.phat(4),Pinj_dist.phat(5));
    hold on;
    % plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
    %     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
    plot(fitx,fitv,'linewidth',2)
    legend('Histogram','Mixture Normal')
    xlabel('Pinj Distance to Source [hops/max hop]')
    ylabel('Density')
    set(gca,'FontSize',20)
else
    [v,edges] = histcounts(pinj_hopfrac,'normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end),max(length(v),100));
    fitv = Pinj_dist.norm_bi(fitx,Pinj_dist.phat(1),Pinj_dist.phat(2),...
            Pinj_dist.phat(3),Pinj_dist.phat(4),Pinj_dist.phat(5));
end

%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
KL = kld(v,Pinj_dist.norm_bi(centers,Pinj_dist.phat(1),Pinj_dist.phat(2),...
            Pinj_dist.phat(3),Pinj_dist.phat(4),Pinj_dist.phat(5)),centers);
if show_fig
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f',KL);
    text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end