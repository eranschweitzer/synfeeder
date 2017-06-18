function KL = kld(p,q,klx)
% The KL divergence D_KL(p||q) = integral_{inf}^{inf} p(x)ln(p(x)/q(x)) dx
% INPUTS
%   v: Emperical function values [p(x)]. Must be a vector
%   klx: p(klx) = v. The spacing of klx is assumed to be unifrom
if length(size(p)) ~=2
    error('Input p must be a vector!')
end
if length(size(p)) ~=2
    error('Input q must be a vector!')
end
if length(size(klx)) ~=2
    error('input klx must be a vector!')
end
if size(p,1) ~= 1
   p = p.'; %ensure row vector 
end
if size(q,1) ~= 1
   q = q.'; %ensure row vector 
end
if size(klx,1) ~= 1
    klx = klx.';
end
v0_mask = p~=0; %need to remove null values
dx = klx(2)-klx(1);
KL = p(v0_mask)*log(p(v0_mask)./q(v0_mask)).'.*dx;