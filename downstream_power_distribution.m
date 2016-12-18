function [KL,v,centers,edges,fitx,fitv] = downstream_power_distribution(n,dP_dist,varargin)
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end

if ~any(strcmp(fieldnames(n),'fid'))
    dP = n.pdownstream(n.d_hop >2 & n.pdownstream >0)./sum(n.p(n.p > 0));
else
    dP = [];
    for fid = unique(n.fid).'
        fid_mask = n.fid == fid;
        pos_dp = n.pdownstream >0;
        hop2 = n.d_hop >2;
        p_pos = n.p > 0;
        norm_factor = sum(n.p(fid_mask & p_pos));
        dP = [dP; n.pdownstream(fid_mask & pos_dp & hop2)./norm_factor]; %#ok<AGROW>
    end
end
[v,edges] = histcounts(dP,'normalization','pdf');
centers = edges(1:end-1) + 0.5*diff(edges);
fitx = linspace(edges(1),edges(end),max(length(centers),100));
fitv = pdf(dP_dist,fitx);
KL = kld(v,pdf(dP_dist,centers),centers);

if show_fig
    loglog(centers,v,'o','MarkerFace','b')
    hold on;
    plot(fitx,pdf(dP_dist,fitx),'linewidth',2)
    xlabel('Downstream P_i/P_{total}, i>2')
    ylabel('Density')
    set(gca,'FontSize',20)
    ax = gca;
    str = sprintf('KL-Divergence: %0.4f\n Target max: 0.25',KL);
    text(0.1*ax.XLim(1)+0.9*ax.XLim(2),0.1*ax.YLim(1)+0.9*ax.YLim(2),str,'Fontsize',16,...
        'VerticalAlignment','top','HorizontalAlignment','right')
end
