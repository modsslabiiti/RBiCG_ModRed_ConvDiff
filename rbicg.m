function [x,xt,flag,i,resvec,resvec_tilde,U1,U1_tilde,U0,U0_tilde] = ... 
    rbicg(A,b,c,tol,maxit,x0,x0_tilde,U,U_tilde,build_space,s,k_ideal,M1,M2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recycled biconjugate gradient for real systems.
% flag = 0 implies converged
%      = 1 maximum iteration reached
%      = 2 break down i.e. division by zero
%
% k_ideal = number of Lanczos vectors we want to recycle.
% k       = num of columns of U.
% k1      = num of columns of U1 (related to building the recycled space).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    [C_hat C_tilde_hat] = binormalize(C, C_tilde, k);

    cum_proj = zeros(k,1);
    cum_projt = zeros(k,1);
    
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
rt = c - applyPrecond(transpose_A, xt, M2', M1');
norm_rt = norm(rt);
if abs((rt'*r)/(norm_r*norm_rt)) < 1e-12
    rt = rand(n,1);
end 
rho = rt' * r;

% Set up for the finding the solution of linear system.
flag = 1;
tolb = tol * norm(b); 
tolc = tol * norm(c); 
resvec = zeros(maxit+1,1);         
resvec_tilde = zeros(maxit+1,1);   
resvec(1) = norm_r;    
resvec_tilde(1) = norm_rt;                 

if (build_space == 1)
    % Initializations related to building the recycle space.
    % - V's and V_tilde's (max columns in each is "s").
    % - T1_extra's and T1_tilde_extra's (the first time number of rows would be
    %                                    s+1, and all other times s+2).
    % - B1 and B1_tilde.
    % - Some flags related to the above 3.
    % 
    % Note: For T's and B's as using "qz", so matrix has to be full; they are
    %       supposed to be very small as k and s are small.

    % V's
    v_old = r / norm_r;
    v_tilde_old = rt / (rt'*v_old);     % This is a different scaling.
    V1 = zeros(n,s);
    V1_tilde = zeros(n,s);

    % T's
    T1_extra = zeros(s+2, s);
    T1_tilde_extra = zeros(s+2, s);
    last_T1_extra_element = 0;
    alpha = -1;             % This one is -1 as it is used in denominator.
    betaa = 0;                 
    r_old = r;

    % B's
    if (~isempty(U))
        B1 = zeros(k,s);
        B1_tilde = zeros(k,s);
        proj = 0;
        projt = 0;
        betaat = 0;
    end

    % flags
    index = 1;        % Used to build V's and T's.
    cycle = 1;        % Length of cycle = Lanczos vectors kept.
end

% loop over maxit iterations (unless convergence or failure)
for i = 1 : maxit
    
    if (build_space == 1)
        r_old_old = r_old;
        r_old = r;
        alpha_old = alpha;
        betaa_old = betaa;
        
        if (~isempty(U))        
            proj_old = proj;
            projt_old = projt;
        end
    end
    
    if i == 1
        p = r;
        pt = rt;
    else
        p = r + betaa * p;
        pt = rt + betaat * pt;
    end

    q = applyPrecond(A, p, M1, M2);
    qt = applyPrecond(transpose_A, pt, M2', M1');

    if (~isempty(U))
        proj = C_hat'*q;
        q = q - C * proj;
        projt = C_tilde_hat'*qt;
        qt = qt - C_tilde * projt;
        
        if (build_space == 1)
            B1(:,index) = (proj - betaa*proj_old)/norm_r; 
            B1_tilde(:,index) = (projt - betaat*projt_old)/(rt'*v_old); 
        end
    end
    
    alpha = rho / (pt' * q);
    alphat = conj(alpha);

    x = x + alpha * p;
    xt = xt + alphat * pt;
    
    if (~isempty(U))
        cum_proj = cum_proj + alpha * proj;
        cum_projt = cum_projt + alphat * projt;
    end

    r = r - alpha * q;
    norm_r = norm(r);
    rt = rt - alphat * qt;
    norm_rt = norm(rt);

    rho_old = rho;
    rho = rt' * r;
    betaa = rho / rho_old;
    betaat = conj(betaa);
    if betaa == 0 || isinf(betaa) || betaat == 0 || isinf(betaat)
        flag = 2;
        break
    end
    
    % check for convergence
    resvec(i+1) = norm_r;
    resvec_tilde(i+1) = norm_rt;
    if norm_r <= tolb && (isempty(nonzeros(c)) || norm_rt <= tolc)
        flag = 0;
        break
    end    

    if (build_space == 1)
        % Update the V1 and V1_tilde matrices.
        V1(:, index) = v_old;                     
        V1_tilde(:, index) = v_tilde_old;

        % Obtain the new v and v_tilde.
        v = r/norm_r;
        v_tilde = rt / (rt'*v);          %This is a different normalization.

        % Updates for next iteration corresponding to V1.
        v_old = v;
        v_tilde_old = v_tilde;

        % Previous tridiagonal
        T1_extra_prev = T1_extra; 
        T1_tilde_extra_prev = T1_tilde_extra;

        % Building T1 and T1_tilde matrices.
        [eta, delta, sai] = triColumn (alpha, alpha_old, betaa_old, r, r_old, r_old_old);
        T1_extra(index, index) = eta;                       
        T1_extra(index+1, index) = delta;                   
        T1_extra(index+2, index) = sai;                     

        % Updates for next iteration corresponding to V1 and T_extra's. 
        index = index + 1;

        % Build the recycle space at end of each cycle
        if (mod(i,s) == 0)
            T1_tilde_extra(1,1) =  conj(last_T1_extra_element);
            tLen = length(T1_extra(2:end-1,:));   % This should be a square matrix
            T1_tilde_extra(2:tLen+1,1:tLen) =  T1_extra(2:end-1,:)';
            T1_tilde_extra(tLen+2,tLen) = conj(-(norm(r_old)*betaa)/(norm(r)*alpha));          
            if (isempty(U))
                % First time it is Ritz Vectors
                if (cycle == 1)                
                    [lambda P P_res P_tilde P_tilde_res] = getGenEigenvecs(k_ideal, T1_extra(2:end-1,:), eye(s,s));

                    % Now build the U0 and U0_tilde from this.
                    U0 = V1*P;
                    U0_tilde = V1_tilde*P_tilde;

%                     % Checking the accuracy of the eigenvectors
%                     % i.e. A x - lambda x = ?
%                     % Note: lambda = conj(lambda_tilde)
%                     U0_complex = V1*P_res;
%                     U0_tilde_complex = V1_tilde*P_tilde_res;
%                     res_norm = norm(applyPrecond(A, U0_complex, M1, M2) - U0_complex*lambda, 'fro');
%                     res_tilde_norm = norm(applyPrecond(transpose_A, U0_tilde_complex, M2', M1') - U0_tilde_complex*lambda, 'fro');
%                     fprintf(1, 'Main space residual for ritz vector computation = %14.4e \n', res_norm);
%                     fprintf(1, 'Dual space residual for ritz vector computation = %14.4e \n', res_tilde_norm);
                else
                    V1_extra = [v_extra_start V1 v];
                    V1_tilde_extra = [v_tilde_extra_start V1_tilde v_tilde];

                    H1 = [eye(k1,k1) zeros(k1,s); zeros(s+2,k1) T1_extra];
                    H2 = [eye(k1,k1) zeros(k1,s); zeros(s+2,k1) T1_tilde_extra];
                    
                    W1 = [U1 V1];
                    W1_tilde = [U1_tilde V1_tilde];
                    
                    W1_hat = [C1 V1_extra];
                    W1_tilde_hat = [C1_tilde V1_tilde_extra];
                    WW1 = W1_tilde_hat' * W1_hat;
                    WW2 = W1_tilde_hat' * W1;                    
                    
                    [U0 U0_tilde P P_tilde] = prepareProblem(W1, W1_tilde, WW1, WW2, H1, H2, k_ideal, A, transpose_A, M1, M2);
                end
            else
                % Build the W, V, T and H matrices.
                if (cycle == 1)
                    V1_extra = [V1 v];
                    V1_tilde_extra = [V1_tilde v_tilde];

                    % Remove the first row as there is no corresp.v vector yet.
                    T1_extra(1, :) = [];
                    T1_tilde_extra(1, :) = [];

                    % First time number of rows in T1_extra is s+1.
                    H1 = [eye(k,k) B1; zeros(s+1,k) T1_extra];
                    H2 = [eye(k,k) B1_tilde; zeros(s+1,k) T1_tilde_extra];

                    W1 = [U V1];
                    W1_tilde = [U_tilde V1_tilde];

                    WW1_term_1 = C_tilde'*C;
                    WW1_term_2 = zeros(k, k);           % For next cycle       
                    WW1_term_3 = zeros(k, s+1);
                    WW1_term_4 = zeros(k, k);           % For next cycle  
                    WW1_term_7 = zeros(s+1, k);
                    WW1_term_9 = V1_tilde_extra'*V1_extra;

                    WW2_term_1 = C_tilde'*U;
                    WW2_term_2 = zeros(k, s);
                    WW2_term_3 = zeros(k, k);           % For next cycle
                    WW2_term_4 = zeros(k, s);           % For next cycle      
                    WW2_term_5 = V1_tilde_extra'*U;
                    WW2_term_6 = V1_tilde_extra'*V1;

                    % Used in WW2_term_3
                    WW1_term_8_prev = zeros(s, k);      % For next cycle  
                    Hdim = s+1;                         % For next cycle  
                    k1_old = k;                         % For next cycle

                    WW1 = [WW1_term_1 WW1_term_3; ...
                            WW1_term_7 WW1_term_9];

                    WW2 = [WW2_term_1 WW2_term_2; ...
                            WW2_term_5 WW2_term_6];                
                else
                    V1_extra = [v_extra_start V1 v];
                    V1_tilde_extra = [v_tilde_extra_start V1_tilde v_tilde];

                    % Subsequent times number of rows in T1_extra is s+2.
                    H1 = [zeros(k,k1) B1; eye(k1,k1) zeros(k1,s); zeros(s+2,k1) T1_extra];
                    H2 = [zeros(k,k1) B1_tilde; eye(k1,k1) zeros(k1,s); zeros(s+2,k1) T1_tilde_extra];

                    W1 = [U1 V1];
                    W1_tilde = [U1_tilde V1_tilde];

                    % This is used almost in all terms
                    PN = P*N;
                    P_tildeM = P_tilde*M;

                    % This last term of WW1 done first because it is used by
                    % terms like:WW1_term_6, WW1_term_8, WW2_term_4,
                    % WW2_term_6. It is just a diagonal matrix but need to be
                    % optimized to avoid "n".
                    WW1_term_9 = V1_tilde_extra'*V1_extra;      

                    % WW1_term_1 remains the same always
                    WW1_term_2 = [WW1_term_2 WW1_term_1*B1]*PN;
                    WW1_term_3 = zeros(k, s+2);
                    WW1_term_4 = P_tildeM'*[WW1_term_4; B1_tilde'*WW1_term_1];
                    WW1_term_5 = Sigma;

                    % For WW1_term_6 = C1_tilde'*V1_extra;
                    V1_prod = zeros(Hdim, s+2);
                    V1_prod(Hdim-1,1) = WW1_term_9(1,1);
                    V1_prod(Hdim,2) = WW1_term_9(2,2);
                    WW1_term_6 = P_tildeM'*[zeros(k1_old, s+2); T1_tilde_extra_prev'*V1_prod];

                    WW1_term_7 = zeros(s+2, k);

                    % For WW1_term_8 = V1_tilde_extra'*C1;
                    WW1_term_8 = [zeros(s+2, k1_old) V1_prod'*T1_extra_prev]*PN;     

                    WW2_term_1 = [WW2_term_1 zeros(k, s)]*PN;
                    % WW2_term_2 remains the same always
                    WW2_term_3 = P_tildeM'*[WW2_term_3 WW2_term_4; ...
                                            WW1_term_8_prev T1_tilde_extra_prev'*WW2_term_6]*PN;

                    % For WW2_term_4 = C1_tilde'*V1;                 
                    WW2_term_4 = WW1_term_6(:,2:end-1);

                    WW2_term_5 = V1_tilde_extra'*U1;            % This is anyway's small

                    % For WW2_term_6 = V1_tilde_extra'*V1;
                    WW2_term_6 = WW1_term_9(:,2:end-1);

                    WW1_term_8_prev = WW1_term_8(2:end-1,:);    % Used in WW2_term_3     
                    Hdim = s+2;                                 
                    k1_old = k1;

                    WW1 = [WW1_term_1 WW1_term_2 WW1_term_3; ...
                            WW1_term_4 WW1_term_5 WW1_term_6; ...
                            WW1_term_7 WW1_term_8 WW1_term_9];

                    WW2 = [WW2_term_1 WW2_term_2; ...
                            WW2_term_3 WW2_term_4; ...
                            WW2_term_5 WW2_term_6];                
                end
                [U0 U0_tilde P P_tilde] = prepareProblem(W1, W1_tilde, WW1, WW2, H1, H2, k_ideal, A, transpose_A, M1, M2);
            end
            cycle = cycle + 1;

            C0 = applyPrecond(A, U0, M1, M2);
            C0_tilde = applyPrecond(transpose_A, U0_tilde, M2', M1');        

            % Orthogonalize the C0 and C0_tilde.
            [U1 C1 U1_tilde C1_tilde k1 M Sigma N] = orthogonalize(U0, C0, U0_tilde, C0_tilde);

            % Updates related to V's and T's.
            v_extra_start = V1(:, end);
            v_tilde_extra_start = V1_tilde(:, end);
            last_T1_extra_element = T1_extra(end,end);
            index = 1;
        end
    end
end                                % for i = 1 : maxit

disp('Done with rbicg');

% Truncate the zeros from resvec
resvec = resvec(1:i+1);
resvec_tilde = resvec_tilde(1:i+1);

% Optimization i.e. update the solution at end.
if (~isempty(U))
    x = x - U * cum_proj;
    xt = xt - U_tilde * cum_projt;
end

% Both U1 and U1_tilde would be missing together if no recycle space is
% built. This can happen in number of cases like when the algorithm
% converges in iterations less that s i.e. length of a cycle.
if (build_space == 1 && ~exist('U1', 'var'))
    U0 = [];
    U1 = [];
    U0_tilde = [];
    U1_tilde = [];
end

%% Compute the next column of the tridiagonal matrix
function [eta, delta, sai] = triColumn (alpha, alpha_old, betaa_old, r, r_old, r_old_old)

% Note that at start:
%   1. betaa_old = 0 => eta = 0 and delta = alpha_inv.
%   2. alpha_old = -1 => no divion by zero.
%   3. If alpha = 0 then it is a break down situation.
%   3. r_old ~= 0 (as it is initial residual) => no divison by zero.
%
% At subsequent steps:
%   1. If alpha_old = 0 or alpha = 0 then it is a break down situation.
%   2. r_old ~= 0 at subsequent steps because otherwise the solution has
%   been found.

alpha_inv = 1/alpha;
constants_ratio = betaa_old/alpha_old;

eta = - norm(r_old_old) / norm (r_old) * constants_ratio;
delta = alpha_inv + constants_ratio;          
sai = - norm(r) / norm (r_old) * alpha_inv;

%% Prepare the generalized eigenvalue problem
function [U0 U0_tilde P P_tilde] = prepareProblem(W1, W1_tilde, WW1, WW2, H1, H2, k_ideal, A, transpose_A, M1, M2)
% Get generalized eigenvectors with the k smallest (absolute)
% generalized eigenvalues ("supposedly" returns real
% eigenvectors). 
gen_A = H2'*WW1*H1;
gen_B = H2'*WW2;            
[lambda P P_res P_tilde P_tilde_res] = getGenEigenvecs(k_ideal, gen_A, gen_B);

% Now build the U0 and U0_tilde from this.
U0 = W1*P;
U0_tilde = W1_tilde*P_tilde;

% % Checking the accuracy of the eigenvectors
% % i.e. A x - lambda x = ?
% U0_complex = W1*P_res;
% U0_tilde_complex = W1_tilde*P_tilde_res;
% res_norm = norm(applyPrecond(A, U0_complex, M1, M2) - U0_complex*lambda, 'fro');
% res_tilde_norm = norm(applyPrecond(transpose_A, U0_tilde_complex, M2', M1') - U0_tilde_complex*lambda, 'fro');
% fprintf(1, 'Main space residual for generalized eigen value = %14.4e \n', res_norm);
% fprintf(1, 'Dual space residual for generalized eigen value = %14.4e \n', res_tilde_norm);

