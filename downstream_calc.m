function [n,e] = downstream_calc(n,e)
%% copy load values to downstream
% each node has at least its own load downstream
n.pdownstream = n.p;
n.qdownstream = n.q;
%% processing order
%we want to process the nodes from the furthest all the way towards the
%source
[~,I] = sort(n.d_hop,'descend');
%% calculate downstream values
nlist = n.id(I(1:end-1)); %we leave out the source
for k = nlist.'
    % add to the predecessor node the downstream load of the current node
    pred = n.pred(k);
    n.pdownstream(pred) = n.pdownstream(pred) + n.pdownstream(k);
    n.qdownstream(pred) = n.qdownstream(pred) + n.qdownstream(k);
    
    emask = find((e.f == pred) & (e.t == n.id(k)));
    if length(emask) > 1
        error('Incorrect assignment of predecessor and child')
    else
        e.pdownstream(emask) = n.pdownstream(k);
        e.qdownstream(emask) = n.qdownstream(k);
    end
end