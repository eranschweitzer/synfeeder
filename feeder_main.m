%main script to create feeders
clear variables; close all;
%% Input
% N = 195;
% Stotal = 23;    %[MVA]
% Pinj_total = 3; %[MW]
in_vars = readtable('input_2.csv');
%% Generate and test

for fid = 1:size(in_vars,1)
    N = in_vars.n(fid);
    Stotal = in_vars.mva(fid);
    Pinj_total = in_vars.pinj(fid);
    fprintf('Processing feeder %d of %d with %d nodes\n',fid,size(in_vars,1),N)
    [n,e] = single_feeder_gen(N,Stotal,Pinj_total);
    n.fid = fid*ones(length(n.id),1);
    e.fid = fid*ones(length(e.id),1); 
    if ~exist('n_tot','var')
        n_tot = n;
        e_tot = e;
    else
        n_tot = table2struct([struct2table(n_tot);struct2table(n)],'ToScalar',true);
        e_tot = table2struct([struct2table(e_tot);struct2table(e)],'ToScalar',true);
    end
end
[KL,T_hist,T_fit] = feeder_testing(n_tot,e_tot,'fig',0);

save('/sine/tmp/eran/radial_feeder_results/input_2_results.mat','n_tot','e_tot','v-7.3')
%%
writetable(struct2tab_with_pad(KL),'~/Dropbox/ASU/SINE/Data/feeder_analysis/input_2_KL.csv');
writetable(struct2tab_with_pad(T_hist),'~/Dropbox/ASU/SINE/Data/feeder_analysis/input_2_T_hist.csv');
writetable(struct2tab_with_pad(T_fit),'~/Dropbox/ASU/SINE/Data/feeder_analysis/input_2_T_fit.csv');
