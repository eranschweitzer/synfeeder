%main script to create feeders
clear variables; close all;
%% Input
N = 195;
Stotal = 23;    %[MVA]
Pinj_total = 3; %[MW]
%% Generate and test 
[n,e] = single_feeder_gen(N,Stotal,Pinj_total);
[KL,T_hist,T_fit] = feeder_testing(n,e,'oneline_only',1);