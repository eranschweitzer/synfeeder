%Single Feeder Generation
% clear variables; 
close all;
%% Inputs
N = 195;
Stotal = 23;    %[MVA]
Pinj_total = 3; %[MW]
U_hv = 110;     %[kV]
U_mv = 10;      %[kV]
% max_overload = 1.5; %maximum allowable overload for a cable
% max_length = 20;%[km] %maximum allowable cable length
min_length = 1e-3; %[km] %minimum allowable cable length
max_parallel = 5; %maximum allowable number of parallel cables
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

%Distributions for negative (injection) load
%includes:
% - pd_p_neg_frac (Beta): for fraction of power injection nodes
% - norm_bi (bimodal normal function and phat [p,mu1,sigma1,mu2,sigma2]:
%   for the hop distance relative to the maximum hop distance where power
%   injections occure
% - pd_p_neg (normal): p_i/sum(p_inj) - 1/N_inj 
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
%% testing
%%%%%%%%%%%%%%%%%%%%%
G = graph(e.f,e.t,e.inom);%sqrt(e.pdownstream.^2 + e.qdownstream.^2));
%% Edge labels & Weight
elabel = cell(size(G.Edges,1),1);
for k = 1:length(e.t)
    eidx = findedge(G,e.f(k),e.t(k));
%     str = sprintf('%d x %d',e.inom(k),e.num_parallel(k));
    str = sprintf('%0.3f',e.pdownstream(k));
%     str = sprintf('%0.3f',e.length(k));
    elabel{eidx} = str;
    G.Edges.Weight(eidx) = e.pdownstream(k);
end
%% plot Graph
figure; 
% plot(G,'layout','layered','source',1,...
%     'nodelabel',1e3*sqrt(n.pdownstream.^2 + n.qdownstream.^2)./(sqrt(3)*n.unom),...
%     'edgelabel',G.Edges.Weight) 
% plot(G,'layout','layered','source',1,...
%     'nodelabel',n.id,...
%     'edgelabel',G.Edges.Weight) 
% hG = plot(G,'layout','layered','source',1,...
%     'nodelabel',n.p,...
%     'edgelabel',elabel);
% hG = plot(G,'layout','layered','source',1,...
%     'nodelabel',n.p);
% hG = plot(G,'layout','layered','source',1,...
%     'edgelabel',elabel,...
%     'linewidth',2);
hG = plot(G,'layout','layered','source',1);
set(gca,'Xtick',[],'Ytick',[])
%% mark nodes with injections
% highlight(hG,n.id(n.p<0),'NodeColor','g','MarkerSize',10)
highlight(hG,n.id(n.p<0),'NodeColor','g');
%% marke edges with reverse flow in red
% highlight(hG,G.Edges.EndNodes(G.Edges.Weight<0),'EdgeColor','r')
G.Edges.EdgeColor = repmat(hG.EdgeColor,length(e.id),1);
G.Edges.EdgeColor(G.Edges.Weight<0,:) = repmat([1 0 0],sum(G.Edges.Weight<0),1);
hG.EdgeColor = G.Edges.EdgeColor;
%% scale node size based on power
marker_max= 10;
marker_min = 4;
msize = (marker_max-marker_min)*abs(n.p)./max(abs(n.p)) + marker_min;
hG.MarkerSize = msize;
%% set line width equal to weight
w = sort(abs(G.Edges.Weight));
G.Edges.LWidths = 6.75*abs(G.Edges.Weight)/max(w(end-1)) + 0.25;
G.Edges.LWidths(G.Edges.LWidths>7) = 7.2;
hG.LineWidth = G.Edges.LWidths;
%% stop excecution here
% error('Excecution Paused')
T_hist = table(); %table for exporting data to csv and tikz
T_fit = table();
KL = table();
%% hop distance distribution
figure;
% plot(x_d_hop,d_hist.hop.values,'o','MarkerFaceColor','b');
h = histogram(n.d_hop,-0.5:1:max(n.d_hop)+0.5,'Normalization','pdf');
x_d_hop = h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges);
hold on;
% plot(x_d_hop,nbinpdf(x_d_hop,phat(1),phat(2)),'linewidth',2);
plot(x_d_hop,pdf(hop_dist,x_d_hop),'linewidth',2);
xlabel('Distance To Source [hop]')
ylabel('Density')
legend('Histogram','Negative Binomial Fit')
set(gca,'FontSize',20)
set(gca,'Xlim',[0 max(n.d_hop)])

T_hist.hopcenters = x_d_hop.';
T_hist.hopv = h.Values.';
T_fit.hopx = x_d_hop.';
T_fit.hopfit = pdf(hop_dist,x_d_hop).';
%% KL Divergence of hop distance
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = h.Values~=0; %need to remove null values
klx = x_d_hop;
dx = klx(2)-klx(1);
KLdiv_hop = h.Values(v0_mask)*log(h.Values(v0_mask)./...
                                    pdf(hop_dist,klx(v0_mask))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n <90%%: 0.39',KLdiv_hop);
text(ax.XLim(1),ax.YLim(2),str,'Fontsize',16,...
    'VerticalAlignment','top','HorizontalAlignment','left')   
KL.hop = KLdiv_hop;
%% Load Distribution
figure;
p_pos_mask = n.p>0;
p_pos_dev = (n.p(p_pos_mask)./sum(n.p(p_pos_mask))) - 1/N;
h = histogram(p_pos_dev,'normalization','pdf');
hold on;
% plot(h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges),...
%     pdf(Pload_dist,h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)),'linewidth',2)
plot(linspace(h.BinEdges(1),h.BinEdges(end)),...
    pdf(Pload_dist,linspace(h.BinEdges(1),h.BinEdges(end))),'linewidth',2)
legend('Histogram','t location-scale fit')
xlabel('P_i/P_{total} - 1/N')
ylabel('Density')
set(gca,'FontSize',20)

%pad Table or vectors with nan and add data
if size(T_hist,1) < length(h.Values)
   T_hist = [T_hist;array2table(nan(length(h.Values)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.pposcenters = (h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges)).';
   T_hist.pposv = h.Values.';
else
   nan_pad = nan(1,size(T_hist,1)-length(h.Values));
   T_hist.pposcenters = [h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges), nan_pad].';
   T_hist.pposv = [h.Values, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.pposx = linspace(h.BinEdges(1),h.BinEdges(end)).';
    T_fit.pposfit = pdf(Pload_dist,linspace(h.BinEdges(1),h.BinEdges(end))).';
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.pposx = [linspace(h.BinEdges(1),h.BinEdges(end)), nan_pad].';
    T_fit.pposfit = [pdf(Pload_dist,linspace(h.BinEdges(1),h.BinEdges(end))), nan_pad].';
end
%% KL Divergence of load distribution
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = h.Values~=0; %need to remove null values
klx = h.BinEdges(1:end-1) + 0.5*diff(h.BinEdges);
dx = klx(2)-klx(1);
KLdiv_load = h.Values(v0_mask)*log(h.Values(v0_mask)./...
                                    pdf(Pload_dist,klx(v0_mask))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n <90%%: 3.4',KLdiv_load);
text(ax.XLim(2),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','right')   
KL.load = KLdiv_load;
%% degree distribution
Pk = histcounts(G.degree,0.5:1:max(G.degree)+0.5,'normalization','pdf');
gam_bi = @(x,p,a1,b1,a2,b2) p*gampdf(x,a1,b1) + (1-p)*gampdf(x,a2,b2);
figure;
plot(Pk,'o','MarkerFace','b')
hold on;
plot(linspace(1,max(G.degree)).',...
    gam_bi(linspace(1,max(G.degree)).',...
    deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5)),...
    'linewidth',2)
xlabel('Degree')
ylabel('Density')
set(gca,'Yscale','log')
set(gca,'FontSize',20)
ax = gca;

%pad Table or vectors with nan and add data
if size(T_hist,1) < length(Pk)
   T_hist = [T_hist;array2table(nan(length(Pk)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.degcenters = (1:length(Pk)).';
   T_hist.degv = Pk.';
else
   nan_pad = nan(1,size(T_hist,1)-length(Pk));
   T_hist.degcenters = [1:length(Pk), nan_pad].';
   T_hist.degv = [Pk, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.degx = linspace(1,max(G.degree)).';
    T_fit.degfit = gam_bi(linspace(1,max(G.degree)).',...
              deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5));
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.degx = [linspace(1,max(G.degree)), nan_pad].';
    T_fit.degfit = [gam_bi(linspace(1,max(G.degree)).',...
    deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5)); nan_pad.'];
end
%% KL Divergence of degree distribution
v0_mask = Pk~=0; %need to remove null values
x = 1:length(v0_mask);
dx = x(2)-x(1);
KLdiv_degree = Pk(v0_mask)*log(Pk(v0_mask)./...
    gam_bi(x(v0_mask),deg_dist(1),deg_dist(2),deg_dist(3),deg_dist(4),deg_dist(5))).'.*dx;
str = sprintf('KL-Divergence: %0.4f\n Target max: 0.2',KLdiv_degree);
text(0.1*ax.XLim(1)+0.9*ax.XLim(2),0.1*ax.YLim(1)+0.9*ax.YLim(2),str,'Fontsize',16,...
    'VerticalAlignment','top','HorizontalAlignment','right')
KL.deg = KLdiv_degree;
%% KL Divergence of Downstream Power
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
[v,bin_edges] = histcounts(n.pdownstream(n.d_hop >2 & n.pdownstream >0)./sum(n.p),'normalization','pdf');

v0_mask = v~=0; %need to remove null values
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
dx = bin_edges(2)-bin_edges(1);
KLdiv_dp = v(v0_mask)*log(v(v0_mask)./pdf(dP_dist,x(v0_mask))).'.*dx;
KL.dp = KLdiv_dp;
%% Downstream Power Distribution
figure;
loglog(x,v,'o','MarkerFace','b')
hold on;
plot(x,pdf(dP_dist,x),'linewidth',2)
xlabel('Downstream P_i/P_{total}, i>2')
ylabel('Density')
set(gca,'FontSize',20)
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n Target max: 0.25',KLdiv_dp);
text(0.1*ax.XLim(1)+0.9*ax.XLim(2),0.1*ax.YLim(1)+0.9*ax.YLim(2),str,'Fontsize',16,...
    'VerticalAlignment','top','HorizontalAlignment','right')

%pad Table or vectors with nan and add data
if size(T_hist,1) < length(v)
   T_hist = [T_hist;array2table(nan(length(v)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.dpcenters = x.';
   T_hist.dpv = v.';
else
   nan_pad = nan(1,size(T_hist,1)-length(v));
   T_hist.dpcenters = [x, nan_pad].';
   T_hist.dpv = [v, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.dpx = linspace(bin_edges(1),bin_edges(end)).';
    T_fit.dpfit = pdf(dP_dist,linspace(bin_edges(1),bin_edges(end))).';
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.dpx = [linspace(bin_edges(1),bin_edges(end)), nan_pad].';
    T_fit.dpfit = [pdf(dP_dist,linspace(bin_edges(1),bin_edges(end))), nan_pad].';
end
%% ratio of i_est and inom
mask = e.inom>0;
r = e.i_est(mask)./(e.inom(mask).*e.num_parallel(mask));
[v,bin_edges] = histcounts(r,'normalization','pdf');
figure;
plot(bin_edges(1:end-1) + 0.5*diff(bin_edges),v,'o',...
    bin_edges(1:end-1) + 0.5*diff(bin_edges),...
    pdf(iest_inom_dist,bin_edges(1:end-1) + 0.5*diff(bin_edges)),...
    'linewidth',2)
set(gca,'Yscale','log')
xlabel('i_{est}/i_{nom} for cables')
ylabel('Density')
legend('Data','Exponential distribution')
set(gca,'FontSize',20)

%pad Table or vectors with nan and add data
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
if size(T_hist,1) < length(v)
   T_hist = [T_hist;array2table(nan(length(v)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.inomcenters = x.';
   T_hist.inomv = v.';
else
   nan_pad = nan(1,size(T_hist,1)-length(v));
   T_hist.inomcenters = [x, nan_pad].';
   T_hist.inomv = [v, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.inomx = linspace(bin_edges(1),bin_edges(end)).';
    T_fit.inomfit = pdf(iest_inom_dist,linspace(bin_edges(1),bin_edges(end))).';
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.inomx = [linspace(bin_edges(1),bin_edges(end)), nan_pad].';
    T_fit.inomfit = [pdf(iest_inom_dist,linspace(bin_edges(1),bin_edges(end))), nan_pad].';
end
%% KL Divergence of iest/inom
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = v~=0; %need to remove null values
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
dx = bin_edges(2)-bin_edges(1);
KLdiv_iest_inom = v(v0_mask)*log(v(v0_mask)./pdf(iest_inom_dist,x(v0_mask))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n Target max: 0.125',KLdiv_iest_inom);
text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
KL.inom = KLdiv_iest_inom;
%% nominal current of coincident edges
r = n.inom_max./n.inom_min;
r = r(n.degree~=1 & ~isinf(r));

r1_mask = r==1;
%%
figure;
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

%% histogram of cable nominal current
figure;
histogram(e.inom(e.inom>0),'Normalization','probability');
xlabel('nominal current [A]')
ylabel('probability')
%% Voltage Drop* distribution
figure;
vdrop_mask = e.i_est > 0 & e.inom > 0 & (e.funom==e.tunom);
vdrop = e.i_est(vdrop_mask).*sqrt(e.x(vdrop_mask).^2 + e.r(vdrop_mask).^2)./(1e3*e.tunom(vdrop_mask));
[v,bin_edges] = histcounts(vdrop,'normalization','pdf');
loglog(bin_edges(1:end-1) + 0.5*diff(bin_edges),v,'o')
hold on;
plot(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
    pdf(vdrop_dist,linspace(mean(bin_edges(1:2)),bin_edges(end))),...
    'linewidth',2)
xlabel('|z|[\Omega]\times i [A]\times (u_{nom}[V])^{-1}')
ylabel('Density')

%pad Table or vectors with nan and add data
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
if size(T_hist,1) < length(v)
   T_hist = [T_hist;array2table(nan(length(v)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.vdropcenters = x.';
   T_hist.vdropv = v.';
else
   nan_pad = nan(1,size(T_hist,1)-length(v));
   T_hist.vdropcenters = [x, nan_pad].';
   T_hist.vdropv = [v, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.vdropx = linspace(bin_edges(1),bin_edges(end)).';
    T_fit.vdropfit = pdf(vdrop_dist,linspace(bin_edges(1),bin_edges(end))).';
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.vdropx = [linspace(bin_edges(1),bin_edges(end)), nan_pad].';
    T_fit.vdropfit = [pdf(vdrop_dist,linspace(bin_edges(1),bin_edges(end))), nan_pad].';
end
%% KL Divergence between voltage drop and GP fit
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = v~=0; %need to remove null values
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
dx = bin_edges(2)-bin_edges(1);
KLdiv_vdrop = v(v0_mask)*log(v(v0_mask)./...
                            pdf(vdrop_dist,x(v0_mask))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f',KLdiv_vdrop);
text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
KL.vdrop = KLdiv_vdrop;
%% cable length distribution
figure;
[v,bin_edges] = histcounts(e.length(e.length>0),50,'normalization','pdf');
loglog(bin_edges(1:end-1) + 0.5*diff(bin_edges),v,'o')

cauchyMod = @(data,x0,gamma) (atan(x0/gamma) + pi/2)^(-1).*(gamma./((data-x0).^2 + gamma.^2));
% Fit modified cauchy distribution to data
% cauchyHat = mle(e.length(e.length>0),'pdf',cauchyMod,'start',[0 1],'lower',[-inf eps],'upper',[inf inf]);
cauchyHat = [0.4807, 0.3595];
hold on; 
plot(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
    length_dist.cauchyMod(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
    length_dist.cauchyHat(1),length_dist.cauchyHat(2)))
xlabel('length [km]')
ylabel('Density')

%pad Table or vectors with nan and add data
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
if size(T_hist,1) < length(v)
   T_hist = [T_hist;array2table(nan(length(v)-size(T_hist,1),size(T_hist,2)),...
                               'VariableNames',T_hist.Properties.VariableNames)];
   T_hist.lengthcenters = x.';
   T_hist.lengthv = v.';
else
   nan_pad = nan(1,size(T_hist,1)-length(v));
   T_hist.lengthcenters = [x, nan_pad].';
   T_hist.lengthv = [v, nan_pad].';
end

if size(T_fit,1) < 100
    T_fit = [T_fit;array2table(nan(100-size(T_fit,1),size(T_fit,2)),...
                               'VariableNames',T_fit.Properties.VariableNames)];
    T_fit.lengthx = linspace(mean(bin_edges(1:2)),bin_edges(end)).';
    T_fit.lengthfit = length_dist.cauchyMod(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
                      length_dist.cauchyHat(1),length_dist.cauchyHat(2)).';
else
    nan_pad = nan(1,size(T_fit,1)-100);
    T_fit.lengthx = [linspace(mean(bin_edges(1:2)),bin_edges(end)), nan_pad].';
    T_fit.lengthfit = [length_dist.cauchyMod(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
                       length_dist.cauchyHat(1),length_dist.cauchyHat(2)), nan_pad].';
end
%% KL Divergence of length
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = v~=0; %need to remove null values
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
dx = bin_edges(2)-bin_edges(1);
KLdiv_length = v(v0_mask)*log(v(v0_mask)./...
                length_dist.cauchyMod(x(v0_mask),length_dist.cauchyHat(1),length_dist.cauchyHat(2))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n Target max: 0.3',KLdiv_length);
text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
KL.length = KLdiv_length;
%% reactance distribution
figure;
[v,bin_edges] = histcounts(e.x(e.x>0),50,'normalization','pdf');
loglog(bin_edges(1:end-1) + 0.5*diff(bin_edges),v,'o')

cauchyMod = @(data,x0,gamma) (atan(x0/gamma) + pi/2)^(-1).*(gamma./((data-x0).^2 + gamma.^2));
% Fit modified cauchy distribution to data
% cauchyHat = mle(e.x(e.x>0),'pdf',cauchyMod,'start',[0 1],'lower',[-inf eps],'upper',[inf inf]);
cauchyHat = [0.0446,0.0351];
hold on; 
plot(linspace(mean(bin_edges(1:2)),bin_edges(end)),...
    cauchyMod(linspace(mean(bin_edges(1:2)),bin_edges(end)),cauchyHat(1),cauchyHat(2)))
xlabel('Reactance [\Omega]')
ylabel('Density')
%% KL Divergence of reactance
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
v0_mask = v~=0; %need to remove null values
x = bin_edges(1:end-1) + 0.5*diff(bin_edges);
dx = bin_edges(2)-bin_edges(1);
KLdiv_reactance = v(v0_mask)*log(v(v0_mask)./...
                cauchyMod(x(v0_mask),cauchyHat(1),cauchyHat(2))).'.*dx;
ax = gca;
str = sprintf('KL-Divergence: %0.4f\n Target max: ',KLdiv_reactance);
text(ax.XLim(1),ax.YLim(1),str,'Fontsize',16,...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
KL.x = KLdiv_reactance;