function z = mult_randn(n, sigma, mu)
%%% sample n times from a multivariate normal distribution 
%%% with covariance matrix sigma and mean vector mu
%%% basically copied from the rand help documentation example

d = size(sigma,1);
if nargin == 2
	mu = zeros(1,d);
end
mu = ensure_col_vect(mu).';

R = chol(sigma);

z = repmat(mu,n,1) + randn(n, d)*R;
