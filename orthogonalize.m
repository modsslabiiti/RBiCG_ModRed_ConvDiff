function [U1, C1, U1_tilde, C1_tilde, k1, M, Sigma, N] = orthogonalize(U0, C0, U0_tilde, C0_tilde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orthogonalize the C0 and C0_tilde.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before we do the SVD, we need to normalize the C0 and C0_tilde
% (scaling). As doing C0, need to do U0 as well.
vecs_in_C0 = size(C0,2);
for count = 1:vecs_in_C0
    tempC0 = 1/norm(C0(:,count));
    C0(:,count) = C0(:,count)*tempC0;
    U0(:,count) = U0(:,count)*tempC0;
    tempC0_tilde = 1/norm(C0_tilde(:,count));
    C0_tilde(:,count) = C0_tilde(:,count)*tempC0_tilde;
    U0_tilde(:,count) = U0_tilde(:,count)*tempC0_tilde;
end

% Now build the U1 and U1_tilde from this. Doing the SVD i.e. A =
% U.Sigma.V* Calling the U as M, and V as N because U and V are
% already used for search space and Lanczos vectors respectively.
% That is: C0_tilde*.C0 = M.Sigma.N*
[M, Sigma, N] = svd(C0_tilde' * C0);

% Remove the columns from C1 and C1_tilde that create this trouble
% of 0 singular values (also from U1 and U1_tilde). Removing from
% Sigma, N and M as that removes it from U1, U1_tilde, C1, and
% C1_tilde.
sig = diag(Sigma);
findIdx = find(abs(sig) >= 1e-6);
if length(findIdx) < length(sig)
    N = N(:, findIdx);
    M = M(:, findIdx);
end

% Now update the U1, U1_tilde, C1, and C1_tilde.
U1 = U0 * N;
C1 = C0 * N;
U1_tilde = U0_tilde * M;
C1_tilde = C0_tilde * M;
k1 = size(U1,2);
Sigma = Sigma(1:k1, 1:k1);
