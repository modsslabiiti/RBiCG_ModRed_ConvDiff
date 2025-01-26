%% Main function
function irkaRailPlots(testcase)

addpath('../')
addpath('rail_models')

switch testcase
    case 1
        n = 1357;
        s = 40;
        k = 10;
		droptol = 0.1;
    case 2
        n = 5177;
        s = 50;
        k = 10;
		droptol = 0.05;
    case 3
        n = 20209;
        s = 40;
        k = 20;
		droptol = 0.01;
    case 4
        n = 79841;
        s = 50;
        k = 20;
        droptol = 0.005;
end
A = mmread(sprintf('rail_%d_c60.A', n));
B = mmread(sprintf('rail_%d_c60.B', n));
C = mmread(sprintf('rail_%d_c60.C', n));
E = mmread(sprintf('rail_%d_c60.E', n));
b = B(:,2);
c = C(6,:);
c = c';

r = 6;
colorOrig = {'-r';'-b+'};
colorRecycle = {'--r';':b+'};
rtol = 1e-6;            % RBiCG Tol
itol = 1e-6;            % IRKA Tol
rmaxit = 200;           % RBiCG Max Iteration
imaxit = 4;             % IRKA Max Iteration
irka_iter = 1;

sigma = logspace(-5,0.7,r)';        % Any ways sorts ascending.
fprintf(1,'Real sigma = %14.4e \n', real(sigma));
fprintf(1,'\n');
fprintf(1,'Imag sigma = %14.4e \n', imag(sigma));

for count = 1:r
    x0{count} = zeros(n,1);
    x0_tilde{count} = zeros(n,1);    
end

% Since doing recycling only for the smallest two shifts
U{1} =[];
U_tilde{1} = [];
U{2} =[];
U_tilde{2} = [];
while (irka_iter <= imaxit)

    [V W U U_tilde x0 x0_tilde] = computeVW(n, r, A, E, sigma, droptol, ...
        b, c, rtol, rmaxit, s, k, U, U_tilde, x0, x0_tilde, irka_iter, ...
        colorOrig, colorRecycle);
    
    sigma_old = sigma;
    er = W'*(E*V);
    ar = W'*(A*V);
    sigma = -eig(ar,er);
    
    % Reflect negative shifts to the right half plane.
    sigma = abs(real(sigma)) + 1i*imag(sigma);

    sigma = sort(sigma);
    if abs((norm(sigma) - norm(sigma_old)) / norm(sigma)) <= itol
        disp('Converged ------------');
        break;
    end
    fprintf(1,'Real sigma = %14.4e \n', real(sigma));
    fprintf(1,'\n');
    fprintf(1,'Imag sigma = %14.4e \n', imag(sigma));    
    
    irka_iter = irka_iter + 1;
end
disp('All done');               % display the answer.

%% subfunction computeWVRecycle
function [VOrth WOrth U U_tilde x0 x0_tilde] = ...
    computeVW(n, r, A, E, sigma, droptol, b, c, tol, maxit, s, k, U, ...
    U_tilde, x0, x0_tilde, irka_iter, colorOrig, colorRecycle)

figure (irka_iter);
hold on;
title(sprintf('%d primary systems with s = %d, k = %d, and drop tol = %1.3f', ...
    n, s, k, droptol),'fontsize', 12, 'fontweight', 'b');
xlabel('Number of iterations', 'fontsize', 12, 'fontweight', 'b');
ylabel('log_{10}||r||', 'fontsize', 12, 'fontweight', 'b');
hold off;

figure (irka_iter+100);
hold on;
title(sprintf('%d dual systems with s = %d, k = %d, and drop tol = %1.3f', ...
    n, s, k, droptol),'fontsize', 12, 'fontweight', 'b');
xlabel('Number of iterations', 'fontsize', 12, 'fontweight', 'b');
ylabel('log_{10}||r||', 'fontsize', 12, 'fontweight', 'b');
hold off;

V = [];
W = [];
for i=1:r
    Mat_A = sigma(i)*E - A;
    [Lp,Up] = luinc(Mat_A, droptol);
    bPrecond = Lp\b;
    cPrecond = Up'\c;
    
    if (i == 1 || i == 2)
        [~, ~, ~, ~, resvec, resvec_tilde, U{i}, U_tilde{i}] ...
            = rbicg(Mat_A, bPrecond, cPrecond, tol, maxit, x0{i}, ...
            x0_tilde{i}, U{i}, U_tilde{i}, 1, s, k, Lp, Up);
        
        resvec = resvec/norm(bPrecond);
        resvec_tilde = resvec_tilde/norm(cPrecond);
        resvec(end) = resvec(end) + 1e-16;
        resvec_tilde(end) = resvec_tilde(end) + 1e-16;
        figure (irka_iter);
        hold on;
        plot(log10(resvec), char(colorRecycle(i)), 'Linewidth', 2, 'MarkerSize', 6);
        hold off;
        figure (irka_iter+100);
        hold on;
        plot(log10(resvec_tilde), char(colorRecycle(i)), 'Linewidth', 2, 'MarkerSize', 6);
        hold off;
    end
    
    [x, x_tilde, ~, ~, resvec, resvec_tilde] ...
        = rbicg(Mat_A, bPrecond, cPrecond, tol, maxit, x0{i}, ...
        x0_tilde{i}, [], [], 0, 0, 0, Lp, Up);

    resvec = resvec/norm(bPrecond);
    resvec_tilde = resvec_tilde/norm(cPrecond);
    resvec(end) = resvec(end) + 1e-16;
    resvec_tilde(end) = resvec_tilde(end) + 1e-16;
    if (i == 1 || i == 2)
        figure (irka_iter);
        hold on;
        plot(log10(resvec), char(colorOrig(i)), 'Linewidth', 2, 'MarkerSize', 6);
        hold off;
        figure (irka_iter+100);
        hold on;
        plot(log10(resvec_tilde), char(colorOrig(i)), 'Linewidth', 2, 'MarkerSize', 6);
        hold off;
    end
    
    x0{i} = x;
    x0_tilde{i} = x_tilde;
    x = Up\x;
    x_tilde = Lp'\x_tilde;
    V = [V x];
    W = [W x_tilde];
end

% To avoid the problem of A_r "Matrix is close to singular or badly
% scaled.", need W and V in an orthogonal basis.
[VOrth, ~] = qr(V,0);
[WOrth, ~] = qr(W,0);
if (isreal(V) && isreal(W))
    %we are fine
else
    disp('Trouble here! V and W complex');
end

figure (irka_iter);
hold on;
legend('Shift 1: Recycling BiCG', ...
        'Shift 1: BiCG', ...
        'Shift 2: Recycling BiCG', ...
        'Shift 2: BiCG');
hold off;
figure (irka_iter+100);
hold on;
legend('Shift 1: Recycling BiCG', ...
        'Shift 1: BiCG', ...
        'Shift 2: Recycling BiCG', ...
        'Shift 2: BiCG');
hold off;
