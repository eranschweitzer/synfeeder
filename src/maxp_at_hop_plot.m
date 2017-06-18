function maxp_at_hop_plot(n,pmax_f)

xhop = min(n.d_hop):max(n.d_hop);
max_p = zeros(size(xhop));
for k = xhop
    hop_mask = n.d_hop == k;
    max_p(xhop==k) = max(n.p(hop_mask));
end

stem(xhop,max_p)
xlabel('hop distance')
ylabel('Max MVA')
hold on;
plot(xhop,pmax_f(xhop),'linewidth',2)
legend('data','power law fit')
% set(gca,'Yscale','log')
ylim([0, max(max_p)+5])
set(gca,'FontSize',20)