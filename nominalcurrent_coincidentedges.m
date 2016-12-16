function nominalcurrent_coincidentedges(n)
%% nominal current of coincident edges
r = n.inom_max./n.inom_min;
r = r(n.degree~=1 & ~isinf(r));

r1_mask = r==1;

h = histogram(r(~r1_mask),1:0.2:ceil(max(r)),'normalization','cdf');
xlabel('max(i_{nom})/min(i_{nom}) (\neq 1)');
ylabel('Cummulative Density')
x_lim = get(gca,'Xlim');
x_tick = get(gca,'XTick');
set(gca,'FontSize',20,'Xlim',x_lim,'XTick',x_tick)
ylim([0 1]);
hold on;
for k = [2,3,4]
    idx = find(h.BinEdges<=k,1,'last') - 1;
    h.Values(idx);
    plot([h.BinEdges(1) h.BinEdges(end)],[h.Values(idx) h.Values(idx)],'linewidth',2)
    str = sprintf('<%d: %0.2f',k,h.Values(idx));
    text(h.BinEdges(idx),h.Values(idx),str,'Fontsize',16,'VerticalAlignment','bottom')
end