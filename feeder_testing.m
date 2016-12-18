function [KL,T_hist,T_fit] = feeder_testing(n,e,err,varargin)
%% check input options
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
I = find(strcmp(varargin,'oneline_only'));
if isempty(I)
    oneline_only = 0;
else 
    oneline_only = varargin{I+1};
    if oneline_only
        show_fig = 0;
    end
end
I = find(strcmp(varargin,'all_onelines'));
if isempty(I)
    all_onelines = 0;
else 
    all_onelines = varargin{I+1};
end
%% determine number of feeders
if ~any(strcmp(fieldnames(n),'fid'))
    feeder_number = 1;
    n.fid = ones(length(n.id),1);
    e.fid = ones(length(e.id),1);
else
    feeder_number = length(unique(n.fid));
end
%% import data
data_import;
%% Create & plot Graph
if (show_fig && (feeder_number==1 || all_onelines)) || oneline_only 
    for fid = unique(n.fid).'
        nfid_mask = n.fid == fid;
        efid_mask = e.fid == fid;
        ng = struct();
        for f = fieldnames(n).'
            ng.(f{1}) = n.(f{1})(nfid_mask);
        end
        eg = struct();
        for f = fieldnames(e).'
            eg.(f{1}) = e.(f{1})(efid_mask);
        end
        G = graph(eg.f,eg.t,eg.inom);
        figure;
        feeder_plot(ng,eg,'G',G)
    end
end
%% initialize data structures
T_hist = struct(); %table for exporting data to csv and tikz
T_fit = struct();
KL = struct();
%% hop distance distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.hop,T_hist.hopv,T_hist.hopcenters,T_hist.hopedges,...
    T_fit.hopx,T_fit.hopv] = hop_distance_distribution(n,hop_dist,'fig',show_fig);
%% Load Distribution (modified to allow for multifeeder)
if show_fig
    figure;
end
[KL.load,T_hist.pposv,T_hist.pposcenters,T_hist.pposedges,...
    T_fit.pposx,T_fit.pposv] = load_distribution(n,Pload_dist,'fig',show_fig);
%% Maximum Load at hop (ok for multifeeder)
if show_fig
    figure;
    maxp_at_hop_plot(n,pmax_f)
end
%% degree distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.deg,T_hist.degv,T_hist.degcenters,Thist.degedges,...
    T_fit.degx,T_fit.degv] = degree_distribution(n,deg_dist,'fig',show_fig);
%% Downstream Power Distribution (modified to allow multifeeder)
if show_fig
    figure;
end
[KL.dp,T_hist.dpv,T_hist.dpcenters,T_hist.dpedges,...
    T_fit.dpx,T_fit.dpv] = downstream_power_distribution(n,dP_dist,'fig',show_fig);
%% i_est/inom Distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.inom,T_hist.inomv,T_hist.inomcenters,T_hist.inomedges,...
    T_fit.inomx,T_fit.inomv] = iest_inom_distribution(e,iest_inom_dist,'fig',show_fig);
%% nominal current of coincident edges (ok for multifeeder)
if show_fig
    figure; 
    nominalcurrent_coincidentedges(n)
end
%% histogram of cable nominal current (ok for multifeeder)
if show_fig
    figure;
    histogram(e.inom(e.inom>0),'Normalization','probability');
    xlabel('nominal current [A]')
    ylabel('probability')
end
%% Voltage Drop distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.vdrop,T_hist.vdropv,T_hist.vdropcenters,T_hist.vdropedges,...
    T_fit.vdropx,T_fit.vdropv] = voltagedrop_distribution(e,vdrop_dist,'fig',show_fig);
%% cable length distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.length,T_hist.lengthv,T_hist.lengthcenters,T_hist.lengthedges,...
    T_fit.lengthx,T_fit.lengthv] = length_distribution(e,length_dist,'fig',show_fig);
%% reactance distribution (ok for multifeeder)
if show_fig
    figure;
end
[KL.x,T_hist.xv,T_hist.xcenters,T_hist.xedges,...
    T_fit.xx,T_fit.xv] = reactance_distribution(e,'fig',show_fig);
%% plot error
if show_fig
    figure;
    stem(err_tot.fid,err_tot.P)
    xlabel('Feeder Number')
    ylabel('P_{actual} - P_{total}')
    set(gca,'FontSize',20)
    
    figure;
    stem(err_tot.fid,err_tot.S)
    xlabel('Feeder Number')
    ylabel('S_{actual} - S_{total}')
    set(gca,'FontSize',20)
    
    figure;
    stem(err_tot.fid,err_tot.Pinj)
    xlabel('Feeder Number')
    ylabel('Pinj_{actual} - Pinj_{total}')
    set(gca,'FontSize',20)
end