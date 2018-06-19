function mpc = parallel_branch_join(mpc)
%%% checks for parallel edges in the mpc branch matrix and joins them into
%%% one equivalen branch. Intended use is primarily on radial systems so
%%% the radial algorithms can be employed.
%%% INPUT
%%%      mpc    matpower case with possibly parallel branches
%%% OUTPUT
%%%      mpc    copy of input case but with parallel branches joined
%%% written by Eran Schweitzer (eranschweitzer@gmail.com) June 2018.

define_constants;
n1 = min(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS));
n2 = max(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS));

branch = zeros(size(mpc.branch));
checked = false(size(mpc.branch,1));
nb = 1;
for k = 1:length(n1)
    if checked(k)
        continue
    end
    test = (n1 == n1(k)) & (n2 == n2(k));
    checked = checked | test;
    branch(nb,:) = mpc.branch(k,:); 
    if sum(test) > 1
        % parallel edges
        % R = (1/R1 + 1/R2 ...)^-1
        % X = (1/X1 + 1/X2 ...)^-1
        % B = B1 + B2 + ...
        % Rate = sum of ratings
        branch(nb,BR_R) = 1/sum(1./mpc.branch(test,BR_R));
        branch(nb,BR_X) = 1/sum(1./mpc.branch(test,BR_X));
        branch(nb,BR_B) = sum(mpc.branch(test,BR_B));
        branch(nb,[RATE_A,RATE_B,RATE_C]) = sum(mpc.branch(test,[RATE_A,RATE_B,RATE_C]),1);
    end
    if all(checked)
        break
    else
        nb = nb + 1;
    end
end

mpc.branch = branch(1:nb,:);





