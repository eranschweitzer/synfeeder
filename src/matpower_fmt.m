function mpc = matpower_fmt(n,e,freq)
%%% create a matpower structure representation of the feeder
%%% described by node structure, n, and edge structure, e

if nargin < 3
    freq = 50; %[HZ]
end
omega = 2*pi*freq;
baseMVA = 10; %use 10MVA base instead of 100. This is more common in distribution circuits

N = length(n.id);
zv   = zeros(N,1);
onev = ones(N,1);
% bus matrix: BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN
bus = [n.id, [3;ones(N-1,1)], n.p, n.q, zv, zv, onev, onev, zv, n.unom, onev, 1.1*onev, 0.9*onev];

M = length(e.id);
branch = zeros(sum(e.num_parallel),13);
cnt = 1;
for b = 1:M
	for bb = 1:e.num_parallel(b)
		if e.funom(b) == e.tunom(b)
			br = e.r(b)/(e.funom(b).^2/baseMVA); %[pu] Ohm/ Zbase= Vbase^2/Sbase
			bx = e.x(b)/(e.funom(b).^2/baseMVA); %[pu] Ohm/ Zbase= Vbase^2/Sbase
            bc = (e.funom(b).^2/baseMVA)*b_from_c(e.c(b),omega); %[pu] uF -> S * Zbase= Vbase^2/Sbase
			tap = 0;
            rating = sqrt(3)*e.funom(b)*e.inom(b)*1e-3; %[MVA] sqrt(3)*vnom_ll*inom = snom (3 phase)
		else
			% note: for transformers inom actually contatins the transformer base MVA
			br = e.r(b)*(baseMVA/e.inom(b)); % convert to system pu 
			bx = e.x(b)*(baseMVA/e.inom(b)); % convert to system pu
			bc = e.c(b); %this should be zero
			% off-nominal tap is ratio of transformer hv rating and system voltage
			tap = (e.funom(b)/n.unom(e.f(b)))/(e.tunom(b)/n.unom(e.t(b))) ; 
            rating = e.inom(b);
		end

		% branch matrix: F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX
		new_br =        [e.f(b), e.t(b), br, bx,   bc,   rating, 0,     0,       tap, 0,     1,         -360, 360];
		branch(cnt,:) = new_br;
		cnt = cnt + 1;
	end
end
	  %GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN
      %PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF
gen = [1, sum(n.p), sum(n.q), 10*sum(n.q), -10*(sum(n.q)), 1.0, baseMVA, 1, 10*sum(n.p), 0, zeros(1,11)];

mpc = struct('baseMVA', baseMVA, 'bus', bus, 'branch', branch, 'gen', gen);
