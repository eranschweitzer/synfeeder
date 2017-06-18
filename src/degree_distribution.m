function [KL,v,centers,edges,fitx,fitv] = degree_distribution(n,deg_dist,varargin)

I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
edges = 0.5:1:max(n.degree)+0.5;
centers = edges(1:end-1) + 0.5*diff(edges);
v = histcounts(n.degree,edges,'normalization','pdf');
fitx = linspace(1,max(n.degree)).';
gam_bi = @(x,p,a1,b1,a2,b2) p*gampdf(x,a1,b1) + (1-p)*gampdf(x,a2,b2);
fitv = gam_bi(linspace(1,max(n.degree)).',...
        deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5));
if show_fig
    plot(v,'o','MarkerFace','b')
    hold on;
    plot(fitx,fitv,'linewidth',2)
    xlabel('Degree')
    ylabel('Density')
    set(gca,'Yscale','log')
    set(gca,'FontSize',20)
    ax = gca;
end

%% KL Divergence of degree distribution
KL = kld(v,gam_bi(centers,deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5)),centers);
if show_fig
    str = sprintf('KL-Divergence: %0.4f\n Target max: 0.2',KL);
    text(0.1*ax.XLim(1)+0.9*ax.XLim(2),0.1*ax.YLim(1)+0.9*ax.YLim(2),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end