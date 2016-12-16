function [KL,T_hist,T_fit] = feeder_testing(n,e,varargin)
%% check input options
I = find(strcmp(varargin,'fig'));
if isempty(I)
    show_fig = 1;
else 
    show_fig = varargin{I+1};
end
I = find(strcmp(varargin,'oneline_only'));
if isempty(I)
    oneline_only = 1;
else 
    oneline_only = varargin{I+1};
    if oneline_only
        show_fig = 0;
    end
end
%% import data
data_import;
%% Create Graph
G = graph(e.f,e.t,e.inom);%sqrt(e.pdownstream.^2 + e.qdownstream.^2));
%% plot feeder
if show_fig || oneline_only
    feeder_plot(n,e,'G',G)
end
%% initialize data structures
T_hist = struct(); %table for exporting data to csv and tikz
T_fit = struct();
KL = struct();
%% hop distance distribution
if show_fig
    figure;
end
[KL.hop,T_hist.hopv,T_hist.hopcenters,T_hist.hopedges,...
    T_fit.hopx,T_fit.hopv] = hop_distance_distribution(n,hop_dist,'fig',show_fig);
%% Load Distribution
if show_fig
    figure;
end
[KL.load,T_hist.pposv,T_hist.pposcenters,T_hist.pposedges,...
    T_fit.pposx,T_fit.pposv] = load_distribution(n,Pload_dist,'fig',show_fig);
%% degree distribution
if show_fig
    figure;
end
[KL.deg,T_hist.degv,T_hist.degcenters,Thist.degedges,...
    T_fit.degx,T_fit.degv] = degree_distribution(G,deg_dist,'fig',show_fig);
%% Downstream Power Distribution
if show_fig
    figure;
end
[KL.dp,T_hist.dpv,T_hist.dpcenters,T_hist.dpedges,...
    T_fit.dpx,T_fit.dpv] = downstream_power_distribution(n,dP_dist,'fig',show_fig);
%% i_est/inom Distribution
if show_fig
    figure;
end
[KL.inom,T_hist.inomv,T_hist.inomcenters,T_hist.inomedges,...
    T_fit.inomx,T_fit.inomv] = iest_inom_distribution(e,iest_inom_dist,'fig',show_fig);
%% nominal current of coincident edges
if show_fig
    figure; 
    nominalcurrent_coincidentedges(n)
end
%% histogram of cable nominal current
if show_fig
    figure;
    histogram(e.inom(e.inom>0),'Normalization','probability');
    xlabel('nominal current [A]')
    ylabel('probability')
end
%% Voltage Drop* distribution
if show_fig
    figure;
end
[KL.vdrop,T_hist.vdropv,T_hist.vdropcenters,T_hist.vdropedges,...
    T_fit.vdropx,T_fit.vdropv] = voltagedrop_distribution(e,vdrop_dist,'fig',show_fig);
%% cable length distribution
if show_fig
    figure;
end
[KL.length,T_hist.lengthv,T_hist.lengthcenters,T_hist.lengthedges,...
    T_fit.lengthx,T_fit.lengthv] = length_distribution(e,length_dist,'fig',show_fig);
%% reactance distribution
if show_fig
    figure;
end
[KL.x,T_hist.xv,T_hist.xcenters,T_hist.xedges,...
    T_fit.xx,T_fit.xv] = reactance_distribution(e,'fig',show_fig);