function [x,flag,i,resvec] = rcgs(A,b,c,tol,maxit,x0,x0_tilde,U,U_tilde,M1,M2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recycled conjugate gradients squared.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Avoid transposing all the time. Also, the matrix would be square.
transpose_A = A';
n = length(A);

% Related to using the previously passed recycled space.
if (~isempty(U))
    % Build the C matrices.
    C = applyPrecond(A, U, M1, M2);
    C_tilde = applyPrecond(transpose_A, U_tilde, M2', M1');
    
    % Orthogonalize the C and C_tilde (should have been done before, but
    % if not, then do it now).
    [U C U_tilde C_tilde k] = orthogonalize(U, C, U_tilde, C_tilde);

    % Normalize columns of C and C_tilde with c1*c1_tilde etc.
    [C_hat, C_tilde_hat] = binormalize(C, C_tilde, k);
    
    % Update the initial guess (get component from U*C_hat and
    % U_tilde*C_tilde_hat).
    r = b - applyPrecond(A, x0, M1, M2);
    temp_proj = C_hat'*r;
    x0 = x0 + U*temp_proj;
    rt = c - applyPrecond(transpose_A, x0_tilde, M2', M1');
    temp_proj = C_tilde_hat'*rt;
    x0_tilde = x0_tilde + U_tilde*temp_proj;
end

% Inital Residuals. If by chance (r,rt) is equal to 0, then take rt as
% random.
x = x0;
r = b - applyPrecond(A, x, M1, M2);
norm_r = norm(r);
xt = x0_tilde;
rt0 = c - applyPrecond(transpose_A, xt, M2', M1');
norm_rt0 = norm(rt0);
if abs((rt0'*r)/(norm_r*norm_rt0)) < 1e-12
    rt0 = rand(n,1);
end 

% Set up for the finding the solution of linear system.
flag = 1;
tolb = tol * norm(b); 
resvec = zeros(maxit+1,1);         
resvec(1) = norm_r;    

% loop over maxit iterations (unless convergence or failure)
for i = 1 : maxit

    rho = rt0'*r;

    if rho == 0
        flag = 2;
        disp('Method Fails');
        break;
    end

    if i == 1
        u = r;
        p = u;
    else
        betaa = rho/rho_old;
        u = r + betaa*q;
        p = u + betaa*(q + betaa*p);
    end

    v = applyPrecond(A, p, M1, M2);
    if (~isempty(U))
        proj = C_hat'*v;
        v = v - C * proj;
    end
    
    alpha = rho/(rt0'*v);
    q = u - alpha*v;
    u_hat = u + q;
    x = x + alpha*u_hat;

    u_hat1 = applyPrecond(A, u_hat, M1, M2);
    if (~isempty(U))
        proj = C_hat'*u_hat1;
        u_hat1 = u_hat1 - C * proj;
        x = x - alpha*(U*proj);
    end
    r = r - alpha*u_hat1;

    norm_r = norm(r);
    resvec(i+1) = norm_r;
    if norm_r <= tolb
        flag = 0;
        break
    end

    rho_old = rho;
end                                % for i = 1 : maxit

% truncate the zeros from resvec
resvec = resvec(1:i+1);

