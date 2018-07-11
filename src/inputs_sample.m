function [N, Stotal, Pinj_total] = inputs_sample(n, use_pinj)
%%%inputs_sample Generates n input samples to the feeder algorithm.
%%%     s = inputs_sample(n, use_pinj);
%%%     [N, Stotal, Pinj_total] = inputs_sample(n, use_pinj);
%%%
%%% s is an n x 3 array, where n is the specified number of samples
%%% and 3 are the inputs to the algorithms [N Stotal Pinj_total].
%%% if use_pinj is set to false then Pinj_total will be forced to zero
%%% for all samples.

inputkde = loadinputkde();
%% sample
s = zeros(n,3);
ptr = 0;
while true
	tmp = round(inputkde.resample(n));
	N = tmp(:,1); Stotal = tmp(:,2); 
	if use_pinj
		Pinj_total = tmp(:,3);
	else
		Pinj_total = zeros(n,1);
	end
	test = (N > 5) & (Stotal > 0) & (Pinj_total <= 0) & (abs(Stotal) > abs(Pinj_total));
	s(ptr + 1:ptr + sum(test),:) = [N(test), Stotal(test), abs(Pinj_total(test))];
	ptr = ptr + sum(test);
	n   = n   - sum(test);
	if n == 0
		if ptr ~= size(s,1)
			error('inputs_sample: something went wrong, n=0 but ptr is not the size of the desired output')
		end
		break
	end
end	

switch nargout
    case 0
        N = s;
    case 1
        N = s;
    case 3
        N = s(:,1); Stotal = s(:,2); Pinj_total = s(:,3);
    otherwise
        error('inputs_sample: outputs must be 1 or 3. nargout=0 is treated as 1.')
end

