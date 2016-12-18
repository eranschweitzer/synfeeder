% This script loads all the necessary data for the radial feeder generation
% algorithm
%% load data libraries
cable_types = load('cable_types.mat','c');
cable_types = cable_types.c;
%% import distributions
%power factor distribution
pf_dist = load('pf_distribution.mat','pf_dist');
pf_dist = pf_dist.pf_dist;
pf_cdf = [pf_dist(:,1) cumsum(pf_dist(:,2))];
pf_cdf(end,2) = 1; %make sure the last bin is all the way to 100%

%nodes at d_hop distribution (Negative Binomial)
hop_dist = load('d_hop_distribution.mat','pd');
hop_dist = hop_dist.pd;

%Distribution of positive load/total load - 1/N (tLocation-Scale)
Pload_dist = load('Pload_distribution.mat','pd');
Pload_dist = Pload_dist.pd;

%Coefficients for maximum load function at a given hop level
%should be used in a function as f(hop) = pmax_coeff(1)*hop^(pmax_coeff(2))
pmax_coeff = load('pmax_coeff','maxp_dhop');
pmax_coeff = pmax_coeff.maxp_dhop;
pmax_f = @(x) pmax_coeff(1)*x.^(pmax_coeff(2)); %convenience function

%Distributions for negative (injection) load
%includes:
% - pd_p_neg_frac (Beta): for fraction of power injection nodes
% - norm_bi (bimodal normal function and phat [p,mu1,sigma1,mu2,sigma2]:
%   for the hop distance relative to the maximum hop distance where power
%   injections occure
% - pd_p_neg (normal): p_i/sum(p_inj) - 1/N_inj 
% - pd_p_neg_degree (poisson): frequency of injections at a node of degree k
% - p_neg_at_h1: fraction of feeders where the root node has an injection
Pinj_dist = load('Power_injection_distributions.mat');

%Distribution of fraction of no load nodes (Beta)
noloadfrac_dist = load('NoLoadFraction_distribution.mat','pd');
noloadfrac_dist = noloadfrac_dist.pd;

%Bimodal Poisson Distribution of no-load at d_hop distance: [p,mu1,mu2]
noloadhop_dist = load('NoLoadAtHop_Param.mat','poiss_bi_phat');
noloadhop_dist = noloadhop_dist.poiss_bi_phat;

%Bimodal Gamma Distribution for node degree: [p,a1,b1,a2,b2]
deg_dist = load('degree_distribution','phat');
deg_dist = deg_dist.phat;

%Maximum degree power law fit [b,a]: a*x^b
deg_max = load('degree_max','pfit');
deg_max = deg_max.pfit;

%Downstream Power Distribution (Generalized Pareto)
dP_dist = load('pdownstream_distribution.mat','pd');
dP_dist = dP_dist.pd;

%Distribution of i_est/i_nom (exponential)
iest_inom_dist = load('iest_inom_distribution.mat','pd');
iest_inom_dist = iest_inom_dist.pd;

%Maximum i_nom at d_hop [d_hop maxinom]: if edge_{d_hop} >= d_hop then max(edge_{inom}) <= maxinom
maxinom_dhop = load('maxinom_dhop.mat','maxinom_dhop');
maxinom_dhop = maxinom_dhop.maxinom_dhop;

%CDF for coincident nominal currents on a node [max(inom)/min(inom) cdf]
%the last value does not have to be one, where it is understood that the
%remaining part is for all ratios exeeding the last entry
coincident_inom_ratio = load('coincident_inom_ratio_cdf.mat','coincident_inom_ratio');
coincident_inom_ratio = coincident_inom_ratio.coincident_inom_ratio;

%Per unit voltage drop distribution for cables
vdrop_dist = load('vdrop_distribution','pd_est');
vdrop_dist = vdrop_dist.pd_est;

%Cable Length Distribution (modified Cauchy)
length_dist = load('length_distribution.mat','cauchyMod','cauchyHat');
%modified cauchy quartile function
length_dist.Q = @(x,x0,gamma) gamma*tan((x-1).*(atan(x0/gamma) + pi/2) + pi/2) + x0;
%Maximum length exponential fit: a*exp(b*x) where x is the hop distance
length_max = load('length_max.mat','fitobject');
length_max = length_max.fitobject;