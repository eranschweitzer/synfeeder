function feeder_plot(n,e,varargin)
% plot feeder from node structure n and edge structure e

I = find(strcmp(varargin,'G'));
if ~isempty(I)
    G = varargin{I+1};
else
    %weight assignment is somewhat arbitrary it is just to create this entry in
    %the graph
    G = graph(e.f,e.t,e.inom);
end
I = find(strcmp(varargin,'label'));
if ~isempty(I)
    print_labels = varargin{I+1};
else
    print_labels = false;
end
%% Edge labels & Weight
elabel = cell(size(G.Edges,1),1);
if isfield(e, 'phasing')
    lstyle = repmat({'-'}, length(e.id),1);
end
for k = 1:length(e.t)
    eidx = findedge(G,e.f(k),e.t(k));
%     str = sprintf('%d x %d',e.inom(k),e.num_parallel(k));
    str = sprintf('%0.3f',e.pdownstream(k));
%     str = sprintf('%0.3f',e.length(k));
    elabel{eidx} = str;
    G.Edges.Weight(eidx) = e.pdownstream(k);
    
    if isfield(e, 'phasing')
        switch e.phasing{k}
            case 'A'
                lstyle{eidx} = '--';
            case 'B'
                lstyle{eidx} = ':';
            case 'C'
                lstyle{eidx} = '-.';
        end
    end
end
%% plot Graph
% plot(G,'layout','layered','source',1,...
%     'nodelabel',1e3*sqrt(n.pdownstream.^2 + n.qdownstream.^2)./(sqrt(3)*n.unom),...
%     'edgelabel',G.Edges.Weight) 
% plot(G,'layout','layered','source',1,...
%     'nodelabel',n.id,...
%     'edgelabel',G.Edges.Weight) 
% hG = plot(G,'layout','layered','source',1,...
%     'nodelabel',n.p,...
%     'edgelabel',elabel);
% hG = plot(G,'layout','layered','source',1,...
%     'nodelabel',n.p);
% hG = plot(G,'layout','layered','source',1,...
%     'edgelabel',elabel,...
%     'linewidth',2);
hG = plot(G,'layout','layered','source',1);
set(gca,'Xtick',[],'Ytick',[])
%% mark nodes with injections
% highlight(hG,n.id(n.p<0),'NodeColor','g','MarkerSize',10)
highlight(hG,n.id(n.p<0),'NodeColor','g');
%% marke edges with reverse flow in red
% highlight(hG,G.Edges.EndNodes(G.Edges.Weight<0),'EdgeColor','r')
G.Edges.EdgeColor = repmat(hG.EdgeColor,length(e.id),1);
G.Edges.EdgeColor(G.Edges.Weight<0,:) = repmat([1 0 0],sum(G.Edges.Weight<0),1);
hG.EdgeColor = G.Edges.EdgeColor;
%% scale node size based on power
marker_max= 10;
marker_min = 4;
msize = (marker_max-marker_min)*abs(n.p)./max(abs(n.p)) + marker_min;
hG.MarkerSize = msize;
%% set line width equal to weight
w = sort(abs(G.Edges.Weight));
G.Edges.LWidths = 6.75*abs(G.Edges.Weight)/max(w(end-1)) + 0.25;
G.Edges.LWidths(G.Edges.LWidths>7) = 7.2;
hG.LineWidth = G.Edges.LWidths;

%% phasing
if isfield(e, 'phasing')
    hG.LineStyle = lstyle;
end

%% label edges and nodes
if print_labels
    hG.EdgeLabel = elabel;
    hG.NodeLabel = n.p;
end