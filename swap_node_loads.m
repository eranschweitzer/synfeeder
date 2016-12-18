function n = swap_node_loads(n,id1,id2)
%swap the load of the two nodes
ptmp = n.p(id1);
qtmp = n.q(id1);
pftmp = n.pf(id1);
n.p(id1) = n.p(id2);
n.q(id1) = n.q(id2);
n.pf(id1) = n.pf(id2);
n.p(id2) = ptmp;
n.q(id2) = qtmp;
n.pf(id2) = pftmp;