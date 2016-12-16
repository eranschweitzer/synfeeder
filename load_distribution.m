function [KL,v,centers,edges,fitx,fitv] = load_distribution(n,Pload_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

p_pos_mask = n.p>0;
N = length(n.id);
p_pos_dev = (n.p(p_pos_mask)./sum(n.p(p_pos_mask))) - 1/N;
if show_fig
    h = histogram(p_pos_dev,'normalization','pdf');
    v = h.Values;
    edges = h.BinEdges;
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end));
    fitv = pdf(Pload_dist,fitx);
    hold on;
    % plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
    %     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
    plot(fitx,fitv,'linewidth',2)
    legend('Histogram','t location-scale fit')
    xlabel('P_i/P_{total} - 1/N')
    ylabel('Density')
    set(gca,'FontSize',20)
else
    [v,edges] = histcounts(p_pos_dev,'normalization','pdf');
    centers = edges(1:end-1) + 0.5*diff(edges);
    fitx = linspace(edges(1),edges(end));
    fitv = pdf(Pload_dist,fitx);
end

%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
KL = kld(v,pdf(Pload_dist,centers),centers);
if show_fig
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f\n <90%%: 3.4',KL);
    text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
        'VerticalAlignment','bottom','HorizontalAlignment','right')
end
