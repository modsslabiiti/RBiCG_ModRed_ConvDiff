function testConvDiff(precond)

addpath('../')

% In the original Fortran code "ibc(1:4) = 1" implied Dirichlit BC at all 4
% ends i.e. south, west, north, east.

% nx: number of grid points in x-direction (not counting boundary points)
% ny: number of grid points in y-direction (not counting boundary points)
nx = 127;
ny = 127;
n = nx*ny;

% Starting and ending points of the grid (xl, xu, yl, yu).
xl = 0;
xu = 1;
yl = 0;
yu = 1;

% dx: (fixed) mesh width in x-direction; typically 1/(nx+1)
% dy: idem
dx = 1/(1+nx);
dy = 1/(1+ny);

[a, rhs] = disnsy(nx+2, ny+2, xl, xu, yl, yu, dx, dy);

% Truncating because disnsy returns arrays of size (nx+2)*(ny+2) instead of
% nx*ny.
a = a(1:n,:);
rhs = rhs(1:n);

% South = a1, West = a2, East = a4, North = a5.
% sp_matrix = spdiags([a(:,1),a(:,2),a(:,3),a(:,4),a(:,5)],...
%                     [-nx,-1,0,1,nx],...
%                     n,n);
sp_matrix = spdiags([[a(nx+1:end,1);zeros(nx,1)], ...
                    [a(2:end,2);zeros(1,1)], ...
                    a(:,3), ...
                    [zeros(1,1);a(1:(n-1),4)], ...
                    [zeros(nx,1);a(1:(n-nx),5)]], ...
                    [-nx,-1,0,1,nx],...
                    n,n);
size(sp_matrix)
sprank(sp_matrix)

tol = 1e-8;
maxit = 250;
% x0 = rand(n,1);
x0 = ones(n,1)./2;

figure (1);
hold all;
title('Convergence curves','fontsize', 12, 'fontweight', 'b');
xlabel('Number of iterations', 'fontsize', 12, 'fontweight', 'b');
ylabel('log_{10}||r||', 'fontsize', 12, 'fontweight', 'b');

switch precond
    % No preconditioning
    case 0
        C = sp_matrix;
        bt = rhs';
        M1 = [];
        M2 = [];
    % Preconditioning as in BiCGStab paper
    otherwise
        C = sp_matrix;
        % opts.type = 'nofill';
        % opts.milu = 'off';
        opts.type = 'crout';
        opts.droptol = 0.2;
        opts.milu = 'off';
        opts.udiag = 0;
        [Lp,Up] = ilu(sp_matrix,opts);
        % [Lp,Up] = luinc(sp_matrix,'0');
        % options.droptol = 0.06;
        % options.droptol = 0.0001;
        % [Lp,Up] = luinc(sp_matrix,options);
        M1 = Lp;
        M2 = Up;
        bt = Lp\rhs';
end

ct = zeros(n,1);  % c is taken as zero based on van der Vorst page 97.
x0_tilde = x0;
s = 80;
k = 20;

%% Empty recycle space.
[x, xt, flag, iter, resvec, resvec_tilde, U1, U1_tilde, U0, U0_tilde] ...
    = rbicg(C, bt, ct, tol, maxit, x0, x0_tilde, [], [], 1, s, k, M1, M2);
disp(sprintf('BiCG: Zero recycle: flag = %d',flag));
resvec(end) = resvec(end) + 1e-16;
plot(log10(resvec), '-r', 'Linewidth', 2);

% [x, flag, iter, resvec] = rcgs(C, bt, ct, tol, maxit, x0, x0_tilde, [], [], M1, M2);
% disp(sprintf('CGS: Zero recycle: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-b', 'Linewidth', 2);
% 
% [x, flag, iter, resvec] = rbicgstab(C, bt, ct, tol, maxit, x0, x0_tilde, [], [], M1, M2);
% disp(sprintf('BiCGStab: Zero recycle: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-g', 'Linewidth', 2);

%% Now the recycle is built first time and so go again.
U = U1;
U_tilde = U1_tilde;

[x, xt, flag, iter, resvec, resvec_tilde, U1, U1_tilde, U0, U0_tilde] ...
    = rbicg(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, 1, s, k, M1, M2);
disp(sprintf('BiCG: Using the first time recycle space: flag = %d',flag));
resvec(end) = resvec(end) + 1e-16;
plot(log10(resvec), '--bo', 'Linewidth', 2, 'MarkerSize', 6);

% [x, flag, iter, resvec] = rcgs(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, M1, M2);
% disp(sprintf('CGS: Using the first time recycle space: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-ob', 'Linewidth', 2);
%
% [x, flag, iter, resvec] = rbicgstab(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, M1, M2);
% disp(sprintf('BiCGStab: Using the first time recycle space: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-og', 'Linewidth', 2);

%% Now the recycle is built subsequent times and so go again.
U = U1;
U_tilde = U1_tilde;

[x, xt, flag, iter, resvec, resvec_tilde, U1, U1_tilde, U0, U0_tilde] ...
    = rbicg(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, 1, s, k, M1, M2);
disp(sprintf('BiCG: Using subsequent recycle space: flag = %d',flag));
resvec(end) = resvec(end) + 1e-16;
plot(log10(resvec), ':gs', 'Linewidth', 2, 'MarkerSize', 6);

% [x, flag, iter, resvec] = rcgs(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, M1, M2);
% disp(sprintf('CGS: Using subsequent recycle space: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-ob', 'Linewidth', 2);
% 
% [x, flag, iter, resvec] = rbicgstab(C, bt, ct, tol, maxit, x0, x0_tilde, U, U_tilde, M1, M2);
% disp(sprintf('BiCGStab: Using subsequent recycle space: flag = %d',flag));
% resvec(end) = resvec(end) + 1e-16;
% plot(log10(resvec), '-og', 'Linewidth', 2);

%% Legend
% legend(...
%     'BICG: Zero recycle', ...
%     'CGS: Zero recycle', ...
%     'BiCGStab: Zero recycle', ...
%     'BICG: First time', ...
%     'CGS: First time', ...
%     'BiCGStab: First time', ...
%     'BICG: Subsequent', ...
%     'CGS: Subsequent', ...
%     'BiCGStab: Subsequent');

legend(...
    'RBICG: First run', ...
    'RBICG: Second run', ...
    'RBICG: Third run')

% legend(...
%     'CGS', ...
%     'BiCGStab', ...
%     'RCGS', ...
%     'RBiCGStab');

%% Final tasks
if (precond ~= 0)
    x_final = Up\x;
end

