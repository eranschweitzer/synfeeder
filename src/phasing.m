function [n,e] = phasing(n,e)
%%% add phasing to feeder

nl = length(e.id);
nb = length(n.id);
A = sparse([e.f;e.t],[e.t;e.f],1, nb, nb); %adjacency matrix

idx.ABC = false(nl,1);
idx.A   = false(nl,1);
idx.B   = false(nl,1);
idx.C   = false(nl,1);
%% sample fraction of 3-phase lines
ph3 = ph3sample();
ph3 = round(ph3*nl); % convert to integer

%% main corridor
idx.ABC(1) = true; %distribution transformer is 3 phase.

%-- mark any edges with parallel branches as 3 phase
idx.ABC(e.num_parallel > 1) = true;
nodes = unique([e.f(idx.ABC); e.t(idx.ABC)]);
if ~check_connected(A(nodes,nodes))
    %-- connect marked sections
    eids = find(idx.ABC);
    for eid = eids.'
       %-- connect marked selections
       if eid == 1
           continue
       end
       pred = e.f(eid);
       while pred ~= 1
           upstream_e = find((e.t == pred) & (e.f == n.pred(pred)),1);
           idx.ABC(upstream_e) = true;
           pred = e.f(upstream_e);
       end
    end
end

%-- add more edges
nodes = unique([e.f(idx.ABC); e.t(idx.ABC)]);
x = false(nb,1);
x(nodes) = true;
while sum(idx.ABC) < ph3
    chldrn = find(logical(A'*x) - x);    % find all adjacent node to current trunk
    eids = find(ismember(e.t, chldrn));  % find corresponding edge ids
    [~, idxtmp] = max(e.pdownstream(eids)); % find maximum loaded branch
    idx.ABC(eids(idxtmp)) = true;
    x(e.t(eids(idxtmp))) = true;
end

%% assign "rolled" up laterals to phases
chldrn = find(logical(A'*x) - x);    % find all adjacent node to current trunk
p = n.pdownstream(chldrn);
hop = n.d_hop(chldrn);
u = phasing_optimization(p, hop);
x = struct();
for field = {'A','B','C'}
    x.(field{:}) = false(nb,1);
    x.(field{:})(chldrn(u.(field{:}))) = true;

    %% mark all downstream nodes of single phase latteral.
    pred = n.pred(x.(field{:}));
    while true
        xnew = (A'*x.(field{:})) | x.(field{:});
        % remove predecessors of originalting nodes
        xnew(pred) = false; 
        if all(xnew == x.(field{:}))
            break
        else
            x.(field{:}) = xnew;
        end
    end
    
    % mark branches
    idx.(field{:})(ismember(e.t, find(x.(field{:})))) = true;
end

%% test that no branch was doubly assigned
if any((idx.ABC + idx.A + idx.B + idx.C) > 1)
    error('phasing: multiply assigned branches.')
elseif any((idx.ABC + idx.A + idx.B + idx.C) == 0)
    error('phasing: unassigned branches.')
end

%% update structures
e.phasing = repmat({'ABC'}, nl, 1);
n.phasing = repmat({'ABC'}, nb, 1);
for field = {'A', 'B', 'C'}
    e.phasing(idx.(field{:})) = field;
    n.phasing(x.(field{:}))   = field;
end

%% -------- local functions ----------
function f = ph3sample()
%%% returns a fraction of 3 phase lines
f = 0.40 + (0.75-0.40)*rand(); %uniform for now

function bool = check_connected(A)
%%% check whether the subgraph defined by edges in vector v is connected

% BFS
x = false(size(A,1),1);
x(1) = true;
while true
    xnew = (A'*x) | x;
    if all(xnew == x)
        break
    else
        x = xnew;
    end
end

bool = all(x);

function u = phasing_optimization(p,hop)

n = length(p);
total_vars = 0;
vars = struct();
for field = {'ua', 'ub', 'uc'}
    vars.(field{:}) = [total_vars+1,total_vars+n];
    total_vars = total_vars + n;
end
for field = {'pa', 'pb', 'pc', 'tab', 'tac', 'tbc'}
    vars.(field{:}) = total_vars+1;
    total_vars = total_vars + 1;
end
for field = {'ha', 'hb', 'hc', 'thab', 'thac', 'thbc'}
    vars.(field{:}) = total_vars+1;
    total_vars = total_vars + 1;
end

%%% p'u{a,b, c} - P{A,B,C} = 0, 
prob.A = sparse([ones(1,n+1),2*ones(1,n+1), 3*ones(1,n+1)],...
            [vars.ua(1):vars.ua(2), vars.pa,...
             vars.ub(1):vars.ub(2), vars.pb,...
             vars.uc(1):vars.uc(2), vars.pc],...
             [p.', -1, p.', -1, p.',-1],...
             3, total_vars);
prob.rhs = zeros(3,1);
prob.sense = repmat('=',3,1);

%%% hop'u{a,b, c} - hop{A,B,C} = 0, 
prob.A = [prob.A; sparse([ones(1,n+1),2*ones(1,n+1), 3*ones(1,n+1)],...
            [vars.ua(1):vars.ua(2), vars.ha,...
             vars.ub(1):vars.ub(2), vars.hb,...
             vars.uc(1):vars.uc(2), vars.hc],...
             [hop.', -1, hop.', -1, hop.',-1],...
             3, total_vars)];
prob.rhs = [prob.rhs; zeros(3,1)];
prob.sense = vertcat(prob.sense, repmat('=',3,1));

%%% ua + ub + uc = 1
prob.A = [prob.A; sparse([1:n,1:n,1:n], ...
    [vars.ua(1):vars.ua(2), vars.ub(1):vars.ub(2), vars.uc(1):vars.uc(2)],...
    ones(1,3*n), n, total_vars)];
prob.rhs = [prob.rhs; ones(n,1)];
prob.sense = vertcat(prob.sense, repmat('=',n,1));

%%% PA - PB > -tab, PA - PB < tab
%%% PA - PC > -tac, PC - PC < tac
%%% PB - PC > -tbc, PB - PC < tbc
for field = {'ab', 'ac', 'bc'}
prob.A = [prob.A;...
          sparse([1,1,1,2,2,2],...
          repmat([vars.(['p',field{:}(1)]), vars.(['p',field{:}(2)]), vars.(['t',field{:}])],1,2),...
          [-1,1,-1, 1,-1,-1],2,total_vars)];
prob.rhs = [prob.rhs; zeros(2,1)];
prob.sense = vertcat(prob.sense, repmat('<',2,1));
end

%%% HA - HB > -thab, HA - HB < thab
%%% HA - HC > -thac, HC - HC < thac
%%% HB - HC > -thbc, HB - HC < thbc
for field = {'ab', 'ac', 'bc'}
prob.A = [prob.A;...
          sparse([1,1,1,2,2,2],...
          repmat([vars.(['h',field{:}(1)]), vars.(['h',field{:}(2)]), vars.(['th',field{:}])],1,2),...
          [-1,1,-1, 1,-1,-1],2,total_vars)];
prob.rhs = [prob.rhs; zeros(2,1)];
prob.sense = vertcat(prob.sense, repmat('<',2,1));
end

%%% objective = tab + tac + tbc + w*(thab + thac + thbc)
w = sum(p)/sum(hop);
prob.obj = full(sparse([vars.tab, vars.tac, vars.tbc, ...
                        vars.thab, vars.thac, vars.thbc],1,...
                        [ones(3,1), w*ones(3,1)],total_vars,1));

prob.vtype = repmat('c',total_vars,1);
for field = {'ua','ub', 'uc'}
    prob.vtype(vars.(field{:})(1):vars.(field{:})(2)) = 'B';
end

prob.modelsense = 'min';
result = gurobi(prob);

if strcmp(result.status, 'OPTIMAL')
    u.A = find(result.x(vars.ua(1):vars.ua(2)) > 0.5);
    u.B = find(result.x(vars.ub(1):vars.ub(2)) > 0.5);
    u.C = find(result.x(vars.uc(1):vars.uc(2)) > 0.5);
end