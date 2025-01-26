function z = applyPrecond(A, y, L, U)

% L and U would be empty at the same time.
if (isempty(L)) 
    z = A*y;
else
    v = U\y;
    w = A*v;
    z = L\w;
end


