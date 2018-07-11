function [Stotal, Pinj] = ncond_sample(n, kde, delta)
%%% return a sample from kde in N, Stotal, Pinj conditional on number
%%% of nodes N = n. The kde represents a distribution f(N, Stotal, Pinj).
%%% The first step is to create a descrete representation of:
%%% f(Stotal, Pinj | N), done by evalutating the kde over a grid and
%%% normalizing.
%%% f(Stotal, Pinj | N) = f(Pinj | Stotal, N) * f(Stotal | N).
%%% Note that expressed as above the two distributions on the right are
%%% independet and in a single variable and can therefore be sampled
%%% independently using their cumulative distribution.
%%% f(Stotal | N) is obtained by saturating f(Stotal, Pinj | N) over all
%%% Pinj.
%%% Once an Stotal sample is obtained, f(Pinj | Stotal, N) is simply a
%%% column in f(Stotal, Pinj | N), which is again normalized and sampled.
%%%
%%% INPUTS
%%%         n   scalar number of nodes in for the desired sample
%%%         kde kde object to sample
%%%         delta (optional) resolution at which Stotal and Pinj will be
%%%               descretized. Can be either a scalar or a 2x1 vector in
%%%               which case delta(1) is the resolution ofr Stotal and
%%%               delta(2) is for Pinj.
%%% OUTPUTS
%%%         Stotal total MVA sample.
%%%         Pinj   total MW power injection sample.

if nargin < 3
    ds   = 0.1;
    dp   = 0.25;
else
    if isscalar(delta)
        ds = delta;
        dp = delta;
    elseif length(delta) == 2
        ds = delta(1);
        dp = delta(2);
    else
        error('ncond_sample: when supplying argument delta it must be either a scalar or a vector of length 2.')
    end
end

%% setup descretization
sigma = sqrt(diag(kde.covariance));
smax = max(kde.dataset(:,2)) + 3*sigma(2);
pmin = min(kde.dataset(:,3)) - 3*sigma(3);

spts = 0:ds:ceil(smax);
ppts = floor(pmin):dp:0;
[S, P] = meshgrid(spts, ppts); % descretization grid

%% descretization of kde given N = n
fsp_given_n = kde.evaluate([n*ones(numel(S),1), S(:), P(:)]);
fsp_given_n = fsp_given_n/(sum(fsp_given_n(:))*ds*dp); % normalize
fsp_given_n = reshape(fsp_given_n, size(S)); % reshape

%% Stotal Sample given N = n
fs_given_n = sum(fsp_given_n, 1).'*dp; %integrate (ie. sum) along the p direction
fs_given_n = fs_given_n/(ds*sum(fs_given_n)); %normalize

Fs_given_n = cumsum(fs_given_n)*ds; % cumulative distribution

r = rand();
idx = find(r <= Fs_given_n, 1, 'first');
Stotal = spts(idx);

%% Pinj Sample given N=n and Stotal = sampled Stotal
fp_given_sn = fsp_given_n(:, idx);
fp_given_sn = fp_given_sn/(dp*sum(fp_given_sn)); %normalize

Fp_given_sn = cumsum(fp_given_sn)*dp; % cumulative distribution

r = rand();
idx  = find( r <= Fp_given_sn, 1, 'first');
Pinj = ppts(idx);


