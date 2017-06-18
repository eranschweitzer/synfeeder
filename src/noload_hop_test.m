function [KL,v,centers,edges,fitx,fitv] = noload_hop_test(n,noloadhop_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
poiss_bi = @(x,p,lambda1,lambda2) ...
                p*poisspdf(x,lambda1) + (1-p)*poisspdf(x,lambda2);
noload_hop = n.d_hop(n.p==0);
if show_fig
    h = histogram(noload_hop,0.5:1:max(noload_hop)+0.5,'normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = centers;
    fitv = poiss_bi(centers,noloadhop_dist.p,noloadhop_dist.mu1,noloadhop_dist.mu2);
    hold on;
    % plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
    %     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
    plot(fitx,fitv,'-o','linewidth',2)
    legend('Histogram','Mixture Poiss.')
    xlabel('No Load Distance to Source [hops]')
    ylabel('Density')
    set(gca,'FontSize',20)
else
    [v,edges] = histcounts(noload_hop,0.5:1:max(noload_hop)+0.5,'normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = centers;
    fitv = poiss_bi(centers,noloadhop_dist.p,noloadhop_dist.mu1,noloadhop_dist.mu2);
end

%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
KL = kld(v,poiss_bi(centers,noloadhop_dist.p,noloadhop_dist.mu1,noloadhop_dist.mu2),centers);
if show_fig
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f',KL);
    text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end