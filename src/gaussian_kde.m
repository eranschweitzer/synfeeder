classdef gaussian_kde < handle
	% essentially a minimal copy of the scipy implementation of a multivariate gaussian just enabling sampling
	properties
		dataset    %(#data, #dims) ie, observations are row vectors
		bw_method
		d          % number of dimentions
		n          % number of samples in data
		factor     % factor multipying the covariance matrix determined by bw_method
		covariance % covariance matrix of dataset scaled by factor
	end
	methods
		function self = gaussian_kde(dataset, varargin)
			self.bw_method = varargin_parse(varargin,'bw_method', 'scott');

			self.dataset = dataset;
			[self.n, self.d] = size(dataset);
			self.set_bandwidth(self.bw_method)
		end
		function set_bandwidth(self, bw_method)
				if strcmp(bw_method,'scott')
					self.factor = self.n^(-1/(self.d + 4));
				elseif strcmp(bw_method,'silverman')
					self.factor = (self.n*(self.d + 2)/4)^(-1/(self.d + 4));
				elseif isscalar(bw_method)
					self.factor = bw_method;
				else
					error('bw_method should be ''scott'', ''silverman'', or a scalar')
				end
				self.covariance = cov(self.dataset)*self.factor^2;
		end
		function z = resample(self, n)
			sigma = mult_randn(n, self.covariance);
			idx   = randi(self.n, n, 1);
			means = self.dataset(idx,:);
			z = means + sigma;
		end

		function z = single_sample(self, actual_vals, subset)
			if nargin == 1
				actual_vals = false;
			elseif nargin == 2;
				subset = 'all';
			end
			if actual_vals
				if strcmp('all', subset)
					I = randi(self.n);
				else
					I = subset(randi(length(subset)));
				end
				z = self.dataset(I,:);
			else
				z = self.resample(1);
			end
		end
	end
end
