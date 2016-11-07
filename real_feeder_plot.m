%Generate the same type of validation figures as the single feeder
%algorithms for actual database feeders
%% open connection
%%%%%%%%%%%%%%%%%%%

%use this one when at ASU
% conn = database('', 'eranpass', 'jeakieK5Lee7', 'org.postgresql.Driver',...
%     'jdbc:postgresql://db.sine.fulton.asu.edu/eran?ssl=true&sslfactory=org.postgresql.ssl.NonValidatingFactory&');

%use the following when outside of ASU and using port forwarding
conn = database('', 'eranpass', 'jeakieK5Lee7', 'org.postgresql.Driver',...
    'jdbc:postgresql://127.0.0.1/eran?ssl=true&sslfactory=org.postgresql.ssl.NonValidatingFactory&');

setdbprefs('DataReturnFormat','structure') %return values as structure
%% Select feeder
f = 45;%1;
reduced = 'true';
%% get data
sql = ['SELECT row_number() OVER (ORDER BY n_id) idx, n_id, p, q, pdownstream_l pdownstream, qdownstream_l qdownstream ',...
        'FROM feeder_node_analysis WHERE reduced is ' reduced ' and f_id= ' num2str(f)];
n = fetch(conn,sql);
nm = sparse(n.n_id,1,n.idx); %map from n_id to consecutive feeder index 

sql = ['SELECT b_id, f, t, inom, length/1000 length, pdownstream_l pdownstream, qdownstream_l qdownstream, i_est ',...
        'FROM feeder_edge_analysis WHERE reduced is ' reduced ' and f_id= ' num2str(f)];
e = fetch(conn,sql);

sql = ['SELECT source FROM feeder_roots WHERE feeder_number=' num2str(f)];
setdbprefs('DataReturnFormat','numeric') %return values as structure
source = fetch(conn,sql);
setdbprefs('DataReturnFormat','structure') %return values as structure
%% Make Graph
G = graph(full(nm(e.f)),full(nm(e.t)),e.length);
%% Edge labels & Weight
elabel = cell(size(G.Edges,1),1);
for k = 1:length(e.t)
    eidx = findedge(G,full(nm(e.f(k))),full(nm(e.t(k))));
%     str = sprintf('%d x %d',e.inom(k),e.num_parallel(k));
    str = sprintf('%0.3f',e.pdownstream(k));
%     str = sprintf('%0.3f',e.length(k));
    elabel{eidx} = str;
    G.Edges.Weight(eidx) = e.pdownstream(k);
end

%% plot Graph
figure; 
% plot(G,'layout','layered','source',1,...
%     'nodelabel',1e3*sqrt(n.pdownstream.^2 + n.qdownstream.^2)./(sqrt(3)*n.unom),...
%     'edgelabel',G.Edges.Weight) 
% plot(G,'layout','layered','source',1,...
%     'nodelabel',n.id,...
%     'edgelabel',G.Edges.Weight) 
% hG = plot(G,'layout','layered','source',full(nm(source)),...
%     'nodelabel',n.p,...
%     'edgelabel',elabel);
% hG = plot(G,'layout','layered','source',full(nm(source)),...
%     'nodelabel',n.p);
% hG = plot(G,'layout','layered','source',full(nm(source)),...
%     'edgelabel',elabel,...
%     'linewidth',2);
hG = plot(G,'layout','layered','source',full(nm(source)));
set(gca,'Xtick',[],'Ytick',[])

%% mark nodes with injections
highlight(hG,n.idx(n.p<0),'NodeColor','g')

%% marke edges with reverse flow in red
% highlight(hG,G.Edges.EndNodes(G.Edges.Weight<0),'EdgeColor','r')
G.Edges.EdgeColor = repmat(hG.EdgeColor,length(e.b_id),1);
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