function err = load_error_check(n,Ptotal,Stotal,Pinj_total)

P_actual = sum(n.p(n.p>0));
S_actual = sum(sqrt(n.p(n.p>0).^2 + n.q(n.p>0).^2));
Pinj_actual = -sum(n.p(n.p<0));

err = struct('P',(P_actual-Ptotal)/Ptotal,...
            'S',(S_actual-Stotal)/Stotal,...
            'Pinj',(Pinj_actual - Pinj_total)/Pinj_total);
