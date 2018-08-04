function [N, Pinj] = scond_sample(s, kde, dp, nmin)
%%% return a sample from kde in N, Stotal, Pinj conditional on total MVA load
%%% Stotal = s. The kde represents a distribution f(N, Stotal, Pinj).
%%% The first step is to create a descrete representation of
%%% f(N, Pinj | Stotal), done by evalutating the kde over a grid and
%%% normalizing. Then,
%%% f(N, Pinj | Stotal) = f(Pinj | Stotal, N) * f(N | Stotal).
%%% Note that expressed as above the two distributions on the right are
%%% independet and in a single variable and can therefore be sampled
%%% independently using their cumulative distribution.
%%% f(N | Stotal) is obtained by saturating f(N, Pinj | Stotal) over all
%%% Pinj.
%%% Once an N sample is obtained, f(Pinj | Stotal, N) is simply a
%%% column in f(N, Pinj | Stotal), which is again normalized and sampled.
%%%
%%% INPUTS
%%%         s   total feeder load in MVA
%%%         kde kde object to sample
%%%         delta (optional) resolution at which Pinj will be
%%%               descretized, N is always descretized as an integer.
%%%         nmin (optional) minimum number of nodes allowed
%%% OUTPUTS
%%%         N      Number of nodes.
%%%         Pinj   total MW power injection sample.
if nargin < 4
    nmin = 5;
    if nargin < 3
        dp   = 0.25;
    end
end
dn = 1;
%% setup descretization
sigma = sqrt(diag(kde.covariance));
nmax = ceil(max(kde.dataset(:,1)) + 3*sigma(1));
pmin = min(kde.dataset(:,3)) - 3*sigma(3);

npts = 0:dn:nmax;
ppts = floor(pmin):dp:0;
[N, P] = meshgrid(npts, ppts); % descretization grid

%% descretization of kde given N = n
fnp_given_s = kde.evaluate([N(:), s*ones(numel(N),1), P(:)]);
fnp_given_s = fnp_given_s/(sum(fnp_given_s(:))*dn*dp); % normalize
fnp_given_s = reshape(fnp_given_s, size(N)); % reshape

%% Stotal Sample given N = n
fn_given_s = sum(fnp_given_s, 1).'*dp; %integrate (ie. sum) along the p direction
fn_given_s = fn_given_s/(dn*sum(fn_given_s)); %normalize

Fn_given_s = cumsum(fn_given_s)*dn; % cumulative distribution

while true
    r = rand();
    idx = find(r <= Fn_given_s, 1, 'first');
    N = npts(idx);
    if N >= nmin
        break
    end
end

%% Pinj Sample given N=n and Stotal = sampled Stotal
fp_given_sn = fnp_given_s(:, idx);
fp_given_sn = fp_given_sn/(dp*sum(fp_given_sn)); %normalize

Fp_given_sn = cumsum(fp_given_sn)*dp; % cumulative distribution

r = rand();
idx  = find( r <= Fp_given_sn, 1, 'first');
Pinj = ppts(idx);


