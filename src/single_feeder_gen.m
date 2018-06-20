function [n,e,err] = single_feeder_gen(N,Stotal,Pinj_total,varargin)
% INPUTS
%   N: number of nodes
%   Stotal: total MVA consumption
%   Pinj_total: total MW injection

%% check inputs
use_pinj = false;
if nargin < 3
	[N, Stotal, Pinj_total] = inputs_sample(1, use_pinj);
end

I = find(strcmp(varargin,'reverse_mitigate'));
if isempty(I)
    reverse_mitigate = 1;
else 
    reverse_mitigate = varargin{I+1};
end
%% Additional Inputs
U_hv = 110;     %[kV]
U_mv = 10;      %[kV]
% max_overload = 1.5; %maximum allowable overload for a cable
% max_length = 20;%[km] %maximum allowable cable length
min_length = 1e-3; %[km] %minimum allowable cable length
max_parallel = 5; %maximum allowable number of parallel cables
% additive penalty for using parallel trnasformers 
par_xfrm_penalty = 0.1; 
%multiply transformer estimated load by this factor before sizing. 
xfrm_buffer = 1.05; 
%% import data
% loads data used in the function
data_import;
%% Initalize variables
n = node_initialize(N);
%initialize enough edges for a tree
e = edge_initialize(N-1);

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
%% Sort Nodes Based on d_hop
[~,I] = sort(n.d_hop);
for f = fields(n).'
   n.(f{1}) = n.(f{1})(I); 
end
n.id = (1:N).'; %reassign consecutive index
%% Connect Nodes
for k = N:-1:2
    pred = find(n.d_hop == n.d_hop(k) - 1);
    
    %find the predecessor that needs the most edges to reach its assigned
    %degree
    [~,idx] = min(n.degree(pred) - n.degree_assign(pred));
    
    % update node entries
    n.degree(k) = n.degree(k) + 1;
    n.degree(pred(idx)) = n.degree(pred(idx)) + 1;
    n.pred(k) = n.id(pred(idx));
    
    % update edge entries
    e.f(k-1) = pred(idx);
    e.t(k-1) = k;
    e.d_hop(k-1) = n.d_hop(k);
    e.funom(k-1) = n.unom(pred(idx));
    e.tunom(k-1) = n.unom(k);
end
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
if Pinj_total > 0
    Ninj = inf;
    while Ninj >= N-N_noload-1  %ensure there will be at least 1 load node
        Ninj = round(Pinj_dist.pd_p_neg_frac.random*N);
        if Ninj == 0
            Ninj = 1;
        end
    end
%% Mark injection nodes and Assign injection
    hop_max = max(n.d_hop);
    Ninj_offset = 0;
    if rand(1) < Pinj_dist.p_neg_at_h1       
        Ninj_offset = 1;
        %assign an injection to root node
        epsilon = random(Pinj_dist.pd_p_neg);
        while 1/Ninj + epsilon < 0
            % we are only assigning positive injection here (negative load)
            epsilon = random(Pload_dist);
        end
        n.p(2) = -Pinj_total*(1/Ninj + epsilon);
        n.q(2) = n.p(2)*tan(acos(n.pf(2)));
    end
    for k = 1:(Ninj - Ninj_offset)
        % select the node based on its hop distance & degree
        idx = [];
        while isempty(idx)
            % determine the hop distance
            p = rand(1);
            if p <= Pinj_dist.phat(1)
                d_tmp = ceil(hop_max*(Pinj_dist.phat(2) + Pinj_dist.phat(3)*randn(1)));
            else
                d_tmp = ceil(hop_max*(Pinj_dist.phat(4) + Pinj_dist.phat(4)*randn(1)));
            end
            if d_tmp > hop_max
                d_tmp = hop_max;
            end
            % determine the node degree
            deg_tmp = 0;
            while deg_tmp < 1
                deg_tmp = Pinj_dist.pd_p_neg_degree.random;
            end
            %find a node at the desired hop distance
            %and doesn't already have an injection assigned
            idx = find((n.d_hop==d_tmp) & ~n.noload & ~(n.p<0));
            
            if ~isempty(idx)
               %pick the node with the most closely matching degree
               [~,idx2] = min(abs(n.degree(idx) - deg_tmp));
               idx = idx(idx2);
            end
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
       while (1/N + epsilon < 0) || (Ptotal*(1/N + epsilon) > pmax_f(n.d_hop(k)))
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
%% make sure positive cap is not exceeded following normalization
for k = 2:N
    if n.p(k) > 0
        if n.p(k) > pmax_f(n.d_hop(k))
            n.p(k) = pmax_f(n.d_hop(k));
        end
    end
end


%% calculate downstream values
[n,e] = downstream_calc(n,e);
%% split up root node
if n.degree(2) > 20
    nnew = ceil(n.degree(2)/20);% - 1;
    new_nodes = node_initialize(nnew);
    new_edges = edge_initialize(nnew);
    neighbors = e.t(e.f == n.id(2));
    pdownstream_init = n.pdownstream(2);
    for k = 1:nnew
        for f = fieldnames(n).'
            %copy the root node's values
            new_nodes.(f{1})(k) = n.(f{1})(2);
        end
        %make some modifications
        new_nodes.id(k) = max(n.id)+k;
        new_nodes.p(k) = 0;
        new_nodes.q(k) = 0;
        new_nodes.pdownstream(k) = 0;
        new_nodes.qdownstream(k) = 0;
        new_nodes.degree(k) = 1; %connected to root
        new_nodes.d_hop(k) = 2;
        new_nodes.pred(k) = 2;
        
        for f = fieldnames(e).'
            %copy edge connection to source
            new_edges.(f{1})(k) = e.(f{1})((e.f == 1) & (e.t == 2));
        end
        %make some modifications
        new_edges.id(k) = max(e.id) + k;
        new_edges.f(k) = 2;
        new_edges.funom(k) = n.unom(2);
        new_edges.t(k) = new_nodes.id(k);
        new_edges.d_hop(k) = 2;
        new_edges.pdownstream(k) = 0;
        new_edges.qdownstream(k) = 0;
        while (new_nodes.pdownstream(k) < pdownstream_init/(nnew)) && ~isempty(neighbors)
            % add to new node
            new_nodes.pdownstream(k) = new_nodes.pdownstream(k) + ...
                                        n.pdownstream(neighbors(1));
            new_nodes.qdownstream(k) = new_nodes.qdownstream(k) + ...
                                        n.qdownstream(neighbors(1));
%             % remove from old on
%             n.pdownstream(2) = n.pdownstream(2) - n.pdownstream(neighbors(1));
%             n.qdownstream(2) = n.qdownstream(2) - n.qdownstream(neighbors(1));
            % mark change in edge and change predecessor of neighbor
            e.f((e.f==2) & (e.t == neighbors(1))) = new_nodes.id(k);
            n.pred(neighbors(1)) = new_nodes.id(k);
            % increment degree of new node and decrement of old
            new_nodes.degree(k) = new_nodes.degree(k) + 1;
%             n.degree(2) = n.degree(2) - 1;
            %remove this neighbor from list
            neighbors = neighbors(2:end); 
        end
        %update downstream power of new edge
        new_edges.pdownstream(k) = new_nodes.pdownstream(k);
        new_edges.qdownstream(k) = new_nodes.qdownstream(k);
        %degree assign is silly in this case. set it to degree
        new_nodes.degree_assign(k) = new_nodes.degree(k);
    end
%     n.degree_assign(2) = n.degree(2);
    n.degree(2) = nnew + 1;
    n.degree_assign(2) = nnew + 1;
    % adjust hop degrees
    n.d_hop(n.d_hop > 1) = n.d_hop(n.d_hop > 1) + 1;
    e.d_hop(e.d_hop > 1) = e.d_hop(e.d_hop > 1) + 1;
    %update edge between source and original root
%     e.pdownstream((e.f == 1) & (e.t == 2)) = n.pdownstream(2);
%     e.qdownstream((e.f == 1) & (e.t == 2)) = n.qdownstream(2);
    %add new nodes and edges to the normal structures
    for f = fieldnames(n).'
        n.(f{1}) = [n.(f{1}) ; new_nodes.(f{1})];
    end
    for f = fieldnames(e).'
        e.(f{1}) = [e.(f{1}) ; new_edges.(f{1})];
    end
    %recalculate downstream values just for safety
    [n,e] = downstream_calc(n,e);
end
%% reverse current mitigation
if reverse_mitigate
    neg_count = sum(n.pdownstream < 0);
    exit_flag = false;
    neg_offset = 0;
    while neg_count > 0 && ~exit_flag
        n_candidates = n.id((n.p <0) & (n.pdownstream < 0));
        if neg_offset >= length(n_candidates)
            break
        end
        [~,I] = sort(n.pdownstream(n_candidates));
        k = n_candidates(I(1+neg_offset));        
        % find candidates that have greater than or equal degree or at
        % least degree greater than 2 (unless original degree is 2)
        swap_candidates = n.id((n.pdownstream > abs(n.p(k))) & ...
                                (n.degree >= n.degree(k) ) );
        if isempty(swap_candidates)
            swap_candidates = n.id( (n.pdownstream > abs(n.p(k))) & (n.degree > 2) );
        end
        if ~isempty(swap_candidates)
            deg_diff = abs(n.degree(swap_candidates) - n.degree(k));
            hop_diff = abs(n.d_hop(swap_candidates) - n.d_hop(k));
            [~,idx] = sort(deg_diff + hop_diff);
            swap_candidates = swap_candidates(idx);
            for s = swap_candidates.'
                %swap the load of the two nodes
                n = swap_node_loads(n,k,s);
                %recalculate downstream power
                [n,e] = downstream_calc(n,e);
                neg_count_new = sum(n.pdownstream < 0);
                if neg_count_new < neg_count
                    break
                else
                    %swap back
                    n = swap_node_loads(n,k,s);
                    %recalculate downstream power
                    [n,e] = downstream_calc(n,e);
                end
            end
            if neg_count_new < neg_count
                neg_count = neg_count_new;
            else
                neg_offset = neg_offset+1;
                if neg_offset == length(n_candidates)
                    exit_flag = true;
                end
            end
        else
            neg_offset = neg_offset+1;
            if neg_offset == length(n_candidates)
                exit_flag = true;
            end
        end
    end
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
               iest_diff = inom_options - e.i_est(k);
               if r_tmp <= 1
                   %weight by cable distribution 
                   %only for under loaded cables
                   inom_diff = inom_diff./...
                   (cable_types.(['u' num2str(e.funom(k))]).count*ones(1,size(inom_diff,2)));
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
       e.rsamp(k) = r_tmp;
       
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
        e.r(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).r(e.cable_id(k)); %[Ohm]
        e.x(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).x(e.cable_id(k)); %[Ohm]
        e.c(k) = e.length(k)*cable_types.(['u' num2str(e.funom(k))]).c(e.cable_id(k)); %[uF]
    end
end

%% Assign Distribution Transfomer Parameters
% Note that by definition this is the first entry in the edge structure
xfrmsest = sqrt(e.pdownstream(1).^2 + e.qdownstream(1).^2); % estimated apparent power in xfrm
xfrmz    = sqrt(xfrmlib.rpu.^2 + xfrmlib.xpu.^2)*[1,0.5,1/3,0.25];
deltaV   = sqrt(xfrmsest/10*xfrmz); %IMPORTANT! Assuming 10 MVA power base
xfrmopt  = xfrmlib.snom*(1:4);
xfrmtest = xfrmopt - xfrm_buffer*xfrmsest;
xfrmtest(xfrmtest < 0 ) = Inf; % effectively remove all undersized xfrm options
xfrmtest = xfrmtest + par_xfrm_penalty*ones(size(xfrmlib.snom))*(0:3); % add parallel penalty
xfrmtest = xfrmtest + deltaV;
[~,xfrmid] = min(xfrmtest(:));
[e.cable_id(1), e.num_parallel(1)] = ind2sub(size(xfrmtest), xfrmid);

% IMPORTANT: The impedance of the transformer is in *per-unit* not in Ohm like the cables
% similarly, value in inom is NOT a current rating but rather a power rating in MVA
e.r(1) = xfrmlib.rpu(e.cable_id(1)); % [pu]
e.x(1) = xfrmlib.xpu(e.cable_id(1)); % [pu]
e.inom(1) = xfrmlib.snom(e.cable_id(1)); %[MVA]

%% calculate error between inputs and output
err = load_error_check(n,Ptotal,Stotal,Pinj_total);
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
