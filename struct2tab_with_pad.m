function T = struct2tab_with_pad(S)
% convert structure with vectors of unequal length to a table and pad
% excess length with NaNs

l = 0;
for f = fieldnames(S).'
    ltmp = length(S.(f{1}));
    if ltmp > l
        l = ltmp;
    end
    %convert to column vector
    if (size(S.(f{1}),2) > 1) && (size(S.(f{1}),1) == 1)
        % row vector case
        S.(f{1}) = S.(f{1}).';
    elseif (size(S.(f{1}),2) > 1) && (size(S.(f{1}),1) > 1)
        error('All fields must be scalars or vectors')
    end
end

for f = fieldnames(S).'
    pad = l - length(S.(f{1}));
   if pad > 0
      S.(f{1}) = [S.(f{1}); nan(pad,1)];
   end  
end

T = struct2table(S);