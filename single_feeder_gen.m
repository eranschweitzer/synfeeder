function [n,e] = single_feeder_gen(N,Stotal,Pinj_total)
% INPUTS
%   N: number of nodes
%   Stotal: total MVA consumption
%   Pinj_total: total MW injection
% clear variables; close all;
%% Additional Inputs
% N = 195;
% Stotal = 23;    %[MVA]
% Pinj_total = 3; %[MW]
U_hv = 110;     %[kV]
U_mv = 10;      %[kV]
% max_overload = 1.5; %maximum allowable overload for a cable
% max_length = 20;%[km] %maximum allowable cable length
min_length = 1e-3; %[km] %minimum allowable cable length
max_parallel = 5; %maximum allowable number of parallel cables
%% import data
data_import;
%% Initalize variables
n = struct('id',(1:N).','p',zeros(N,1),'q',zeros(N,1),'pf',zeros(N,1),...
            'unom',zeros(N,1),...
            'noload',false(N,1),'d_hop',zeros(N,1),...
            'pdownstream',zeros(N,1),'qdownstream',zeros(N,1),...
            'degree',zeros(N,1), 'degree_assign',zeros(N,1),...
            'inom_min',inf(N,1),'inom_max',zeros(N,1));
%initialize enough edges for a tree
e = struct('id',(1:N-1).','f',zeros(N-1,1),'t',zeros(N-1,1),...
            'funom',zeros(N-1,1),'tunom',zeros(N-1,1),'d_hop',zeros(N-1,1),...
            'pdownstream',zeros(N-1,1),'qdownstream',zeros(N-1,1),...
            'i_est',zeros(N-1,1),'inom',-1*ones(N-1,1),'overload',false(N-1,1),...
            'num_parallel',ones(N-1,1),'cable_id',zeros(N-1,1),...
            'length',zeros(N-1,1),'r',zeros(N-1,1),'x',zeros(N-1,1));

%Ptotal is the total power calculated as the total MVA times the average
%power factor. Qtotal is not calculated since it is not needed for the
%generation
Ptotal = Stotal*sum(prod(pf_dist,2)); 

%% Assign Powerfactor
for k =2:N %skip source since it will not have any load
    p = rand(1);
    idx = find(pf_cdf(:,2) > p,1,'first');
    n.pf(k) = pf_cdf(idx,1);
end

%% Assign hop distances
% node 1 is by default the source at distance 0 (no change)
% node 2 is the root at distance 1
n.d_hop(2) = 1;
for k = 3:N
    hop_tmp = random(hop_dist); %generate a random sample
    while (hop_tmp < 2) && (hop_tmp < N-1)
        % no more 0 or 1 distance nodes allowed
        % no hop distances greater than the total number of nodes
        hop_tmp = random(hop_dist); %generate a random sample
    end
    n.d_hop(k) = hop_tmp;
end
%% Ensure consecutive hop distances
distinct_hop = unique(n.d_hop);
while ~all(diff(distinct_hop)==1);
    hop_threshold = find(diff(distinct_hop)~=1,1,'first');
    n.d_hop(n.d_hop > hop_threshold) = n.d_hop(n.d_hop > hop_threshold) - 1;
    distinct_hop = unique(n.d_hop);
end

%% Assign Nominal Voltage
n.unom(1) = U_hv;
%for now make all nodes the same voltage
n.unom(2:end) = U_mv;
%% Determine Number of No Load Nodes
N_noload = floor(N*random(noloadfrac_dist));
while N_noload < 1
    % at least the source will have no load
    N_noload = floor(N*random(noloadfrac_dist));
end
n.noload(1) = 1; %mark the source as having no load

%% Mark No Load Nodes
for k = 1:N_noload-1 %minus one is for the source node
    idx = [];
    while isempty(idx)
        p = rand(1);
        if p <= noloadhop_dist.p
            d_tmp = poissrnd(noloadhop_dist.mu1);
        else
            d_tmp = poissrnd(noloadhop_dist.mu2);
        end
        idx = find((n.d_hop==d_tmp) & ~n.noload,1);
    end
    n.noload(idx) = 1;
end

%% Determine Number of nodes with injections
if Pinj_total > 0;
    Ninj = ceil(Pinj_dist.pd_p_neg_frac.random*N);
%% Mark injection nodes and Assign injection
    hop_max = max(n.d_hop);
    for k = 1:Ninj
        % select the node based on its hop distance
        idx = [];
        while isempty(idx)
            p = rand(1);
            if p <= Pinj_dist.phat(1)
                d_tmp = ceil(hop_max*(Pinj_dist.phat(2) + Pinj_dist.phat(3)*randn(1)));
            else
                d_tmp = ceil(hop_max*(Pinj_dist.phat(4) + Pinj_dist.phat(4)*randn(1)));
            end
            if d_tmp > hop_max
                d_tmp = hop_max;
            end
            %find a node at the desired hop distance that is not marked as
            %no load and doesn't already have an injection assigned
            idx = find((n.d_hop==d_tmp) & ~n.noload & ~(n.p<0),1);
        end
        
        % assign the power injection 
        epsilon = random(Pinj_dist.pd_p_neg);
        while 1/Ninj + epsilon < 0
            % we are only assigning positive injection here (negative load)
            epsilon = random(Pload_dist);
        end
        n.p(idx) = -Pinj_total*(1/Ninj + epsilon);
        n.q(idx) = n.p(idx)*tan(acos(n.pf(idx)));
    end
end
%% Assign Load
for k = 2:N
    if ~n.noload(k) && ~(n.p(k)<0) %not marked as no load or has a power injection
       epsilon = random(Pload_dist);
       while 1/N + epsilon < 0
           % we are only assigning positive load here
           epsilon = random(Pload_dist);
       end
       n.p(k) = Ptotal*(1/N + epsilon);
       n.q(k) = n.p(k)*tan(acos(n.pf(k)));
    end
end

%% Get Load Closer To Input Value?
%positive load
Pload = sum(n.p(n.p>0));
n.p(n.p>0) = n.p(n.p>0)*(Ptotal/Pload);
n.q(n.p>0) = n.q(n.p>0)*(Ptotal/Pload);
%negative load
Pinj = -1*sum(n.p(n.p<0));
n.p(n.p<0) = n.p(n.p<0)*(Pinj_total/Pinj);
n.q(n.p<0) = n.q(n.p<0)*(Pinj_total/Pinj);
%% Assign degree
%degree is known for source, root and furthest nodes
n.degree_assign(1) = 1;
n.degree_assign(2) = sum(n.d_hop==2) + 1;
n.degree_assign(n.d_hop==max(n.d_hop)) = 1;

for k = 3:N
    if n.degree_assign(k)==0 %degree has not been assigned yet
        deg_tmp = inf;
        while deg_tmp > ceil(deg_max(2)*n.d_hop(k)^deg_max(1))
            p = rand(1);
            if p <= deg_dist(1)
                deg_tmp = random('gamma',deg_dist(2),deg_dist(3));
            else
                deg_tmp = random('gamma',deg_dist(4),deg_dist(5));
            end
            %boosting degree 2 nodes
            if abs(deg_tmp - 2) < 1
                deg_tmp=2;
            else
                deg_tmp = round(deg_tmp);
            end
        end
        n.degree_assign(k) = deg_tmp;
    end
end
%% copy load values to downstream
% each node has at least its own load downstream
n.pdownstream = n.p;
n.qdownstream = n.q;
%% Sort Nodes Based on d_hop
[~,I] = sort(n.d_hop);
for f = fields(n).'
   n.(f{1}) = n.(f{1})(I); 
end
n.id = (1:N).'; %reassign consecutie index
%% Connect Nodes
for k = N:-1:2
    pred = find(n.d_hop == n.d_hop(k) - 1);
    
    %find the predecessor that needs the most edges to reach its assigned
    %degree
    [~,idx] = min(n.degree(pred) - n.degree_assign(pred));
    
    % update node entries
    n.pdownstream(pred(idx)) = n.pdownstream(pred(idx)) + n.pdownstream(k);
    n.qdownstream(pred(idx)) = n.qdownstream(pred(idx)) + n.qdownstream(k);
    n.degree(k) = n.degree(k) + 1;
    n.degree(pred(idx)) = n.degree(pred(idx)) + 1;
    
    % update edge entries
    e.f(k-1) = pred(idx);
    e.t(k-1) = k;
    e.d_hop(k-1) = n.d_hop(k);
    e.funom(k-1) = n.unom(pred(idx));
    e.tunom(k-1) = n.unom(k);
    e.pdownstream(k-1) = n.pdownstream(k);
    e.qdownstream(k-1) = n.qdownstream(k);
end

%% Estimated Current
% i_est [A] = 1e3*s_est[MVA]/(sqrt(3)*unom [kV]) 
% Unom always taken on "to" side
for k = 1:length(e.pdownstream)
    e.i_est(k) = 1e3*sqrt(e.pdownstream(k)^2 + e.qdownstream(k)^2)/... 
                                                      (sqrt(3)*e.tunom(k));
end

%% Determine nominal current for cables with downstream load
for k = length(e.inom):-1:1 %moves from far out towards the source
    if (e.funom(k) == e.tunom(k)) && (e.i_est(k) ~= 0) %not a transformer and some current
       maxinom_tmp = maxinom_dhop(find(e.d_hop(k)>= maxinom_dhop(:,1),1,'last'),2);
       if isempty(maxinom_tmp)
           % no restriction on maximum nominal current (other than
           % available cables of course
           maxinom_tmp = inf;
       end
       
       same_inom_test = rand(1) < (2/3); %try and assign same inom as downstream
       if n.degree(e.t(k)) > 1 && (e.i_est(k) <= n.inom_max(e.t(k))) && same_inom_test
          % assign the same nominal current as maximum downstream branch
          I = find(cable_types.(['u' num2str(e.funom(k))]).inom == n.inom_max(e.t(k)));
       else
           flag = false;
           while ~flag
               par_flag = false;
               while ~par_flag
                   r_tmp = random(iest_inom_dist);  %random sample from the Exp() distribution
                   inom_tmp = e.i_est(k)/r_tmp;
                   if (inom_tmp < max(cable_types.(['u' num2str(e.funom(k))]).inom) || ...
                        e.i_est(k) > max(cable_types.(['u' num2str(e.funom(k))]).inom)) && ...
                        (inom_tmp < max_parallel*max(cable_types.(['u' num2str(e.funom(k))]).inom) || ...
                        e.i_est(k) > max_parallel*max(cable_types.(['u' num2str(e.funom(k))]).inom))
                        % only allow parralel lines when the estimated
                        % current is greater than the maximum available
                        % cable type
                        % only allow the maximum number of parallel lines
                        % unless i_est is greather tha the largest cable in
                        % parallel with itself the maximum number of times
                        par_flag = true;
                   end
               end
               inom_options = cable_types.(['u' num2str(e.funom(k))]).inom;
               par_test = min(floor(inom_tmp./cable_types.(['u' num2str(e.funom(k))]).inom));
               if par_test > 0
                    par_cable = cable_types.(['u' num2str(e.funom(k))]).inom*(2:max(2,par_test));
                    inom_options  = [inom_options par_cable]; %#ok<AGROW>
               end
               
               %coincidence ratio 
%                r_ind = find(rand(1) < coincident_inom_ratio(:,2),1,'first');
%                if ~isempty(r_ind)
%                    max_tf = max([n.inom_max(e.t(k)), n.inom_max(e.f(k))]);
%                    min_tf = min([n.inom_min(e.t(k)), n.inom_min(e.f(k))]);
%                    if max_tf > 0 && ~isinf(min_tf)
%                        coincidence_mask = inom_options <= min_tf*coincident_inom_ratio(r_ind,1);
%                    else
%                        coincidence_mask = true(size(inom_options));
%                    end
%                else
%                    coincidence_mask = true(size(inom_options));
%                end
               
               inom_diff = abs(inom_options - inom_tmp);
               %weight by cable distribution
               inom_diff = inom_diff./...
                   (cable_types.(['u' num2str(e.funom(k))]).count*ones(1,size(inom_diff,2)));
               iest_diff = inom_options - e.i_est(k);
               if r_tmp <= 1
                   inom_diff(iest_diff<0) = inf;
               else
                   inom_diff(iest_diff > 0) = inf;
               end
               %find the closest available nominal current cable
               if ~all(isinf(inom_diff(:)))
                   [~,I] = min(inom_diff(:));
                   [I,npar] = ind2sub(size(inom_diff),I); 
                   inom_tmp = inom_options(I,npar);

                   if (inom_tmp >= n.inom_max(e.t(k))) && ...
                           (inom_options(I,1) < maxinom_tmp)
                       flag = true;
                   end
               end
           end
       end
       %update edge properties
       e.inom(k) = cable_types.(['u' num2str(e.funom(k))]).inom(I);
       e.cable_id(k) = I;
       e.num_parallel(k) = npar;
       e.overload(k) = r_tmp > 1;
       
       %update node properties
       n.inom_min(e.f(k)) = min([n.inom_min(e.f(k)) e.inom(k)]);
       n.inom_min(e.t(k)) = min([n.inom_min(e.t(k)) e.inom(k)]);
       n.inom_max(e.f(k)) = max([n.inom_max(e.f(k)) e.inom(k)]);
       n.inom_max(e.t(k)) = max([n.inom_max(e.t(k)) e.inom(k)]);
    end
end

%% Determine nominal current for cables without downstream load
% Assignment is taken as average of upstream node max and min inom values
idx = find((e.funom == e.tunom) & (e.i_est == 0)); %branch index
for k = idx.'
    inom_tmp = mean([n.inom_max(e.f(k)) n.inom_min(e.f(k))]);
    inom_options = cable_types.(['u' num2str(e.funom(k))]).inom;
    inom_diff = abs(inom_options - inom_tmp);
    iest_diff = inom_options - e.i_est(k);
%     if r_tmp <= 1
%        inom_diff(iest_diff<0) = inf;
%     else
%        inom_diff(iest_diff > 0) = inf;
%     end
    [~,I] = min(inom_diff(:));
    
    %update edge properties
    e.inom(k) = cable_types.(['u' num2str(e.funom(k))]).inom(I);
    e.cable_id(k) = I;
    
    %update node properties (only to side necessary
    n.inom_min(e.t(k)) = min([n.inom_min(e.t(k)) e.inom(k)]);
    n.inom_max(e.t(k)) = max([n.inom_max(e.t(k)) e.inom(k)]);
end

%% Determine Cable Length and Impedance
for k = 1:length(e.inom)
    if (e.funom(k) == e.tunom(k)) && (e.inom(k) > 0) % cable type has been assigned
        lflag = false;
        while ~lflag
            % draw a sample from modified Caucy Distribution
            ltmp = length_dist.Q(rand(1),length_dist.cauchyHat(1),length_dist.cauchyHat(2));
            if (ltmp <= length_max.a*exp(length_max.b*e.d_hop(k))) && ...
                    (ltmp >= min_length)
                % if lenght is within range
                lflag = true;
            end
        end
        e.length(k) = ltmp;
        e.r(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).r(e.cable_id(k));
        e.x(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).x(e.cable_id(k));
    end
end
%% Determine Cable Length and total impedance for cables with i_est>0
% for k = 1:length(e.inom)
%     if (e.funom(k) == e.tunom(k)) && (e.i_est(k) > 0) && (e.inom(k) > 0)
%         lflag = false;
%         while ~lflag
%             vdrop_tmp = random(vdrop_dist);
%             while (e.d_hop(k)>3 && vdrop_tmp > 0.03) || vdrop_tmp > 0.05
%                 vdrop_tmp = random(vdrop_dist);
%             end
%             z = (vdrop_tmp*e.tunom(k)*1e3)/e.i_est(k);
%             z_km = sqrt(...
%                     cable_types.(['u' num2str(e.funom(k))]).r(e.cable_id(k))^2 + ...
%                     cable_types.(['u' num2str(e.funom(k))]).x(e.cable_id(k))^2 ...
%                     );
%             if (z/z_km < max_length) && (z/z_km > min_length)
%                 lflag = true;
%             end
%         end
%         e.length(k) = z/z_km;
%         e.r(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).r(e.cable_id(k));
%         e.x(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).x(e.cable_id(k));
%     end
% end
%% separation between generation and testing code

%%%%%%%%%%%%%%%%%%%%%
% %% testing
% %%%%%%%%%%%%%%%%%%%%%
% G = graph(e.f,e.t,e.inom);%sqrt(e.pdownstream.^2 + e.qdownstream.^2));
% %% plot feeder
% feeder_plot(n,e,'G',G)
% %% stop excecution here
% T_hist = struct(); %table for exporting data to csv and tikz
% T_fit = struct();
% KL = struct();
% %% hop distance distribution
% figure;
% [KL.hop,T_hist.hopv,T_hist.hopcenters,T_hist.hopedges,...
%     T_fit.hopx,T_fit.hopv] = hop_distance_distribution(n,hop_dist,'fig',1);
% %% Load Distribution
% figure;
% [KL.load,T_hist.pposv,T_hist.pposcenters,T_hist.pposedges,...
%     T_fit.pposx,T_fit.pposv] = load_distribution(n,Pload_dist,'fig',1);
% %% degree distribution
% figure;
% [KL.deg,T_hist.degv,T_hist.degcenters,Thist.degedges,...
%     T_fit.degx,T_fit.degv] = degree_distribution(G,deg_dist,'fig',1);
% %% Downstream Power Distribution
% figure;
% [KL.dp,T_hist.dpv,T_hist.dpcenters,T_hist.dpedges,...
%     T_fit.dpx,T_fit.dpv] = downstream_power_distribution(n,dP_dist,'fig',1);
% %% i_est/inom Distribution
% figure;
% [KL.inom,T_hist.inomv,T_hist.inomcenters,T_hist.inomedges,...
%     T_fit.inomx,T_fit.inomv] = iest_inom_distribution(e,iest_inom_dist,'fig',1);
% %% nominal current of coincident edges
% figure; 
% nominalcurrent_coincidentedges(n)
% %% histogram of cable nominal current
% figure;
% histogram(e.inom(e.inom>0),'Normalization','probability');
% xlabel('nominal current [A]')
% ylabel('probability')
% %% Voltage Drop* distribution
% figure;
% [KL.vdrop,T_hist.vdropv,T_hist.vdropcenters,T_hist.vdropedges,...
%     T_fit.vdropx,T_fit.vdropv] = voltagedrop_distribution(e,vdrop_dist,'fig',1);
% %% cable length distribution
% figure;
% [KL.length,T_hist.lengthv,T_hist.lengthcenters,T_hist.lengthedges,...
%     T_fit.lengthx,T_fit.lengthv] = length_distribution(e,length_dist,'fig',1);
% %% reactance distribution
% figure;
% [KL.x,T_hist.xv,T_hist.xcenters,T_hist.xedges,...
%     T_fit.xx,T_fit.xv] = reactance_distribution(e,'fig',1);