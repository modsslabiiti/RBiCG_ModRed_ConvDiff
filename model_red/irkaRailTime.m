%% Main function
function irkaRailTime(testcase, recycleFlag)

addpath('../')
addpath('rail_models')

switch testcase
    case 1
        model_size = 20209;
        s = 40;
        k = 20;
		droptol = 0.01;
    case 2
        model_size = 79841;
        s = 50;
        k = 20;
        droptol = 0.005;
end
A = mmread(sprintf('rail_%d_c60.A', model_size));
B = mmread(sprintf('rail_%d_c60.B', model_size));
C = mmread(sprintf('rail_%d_c60.C', model_size));
E = mmread(sprintf('rail_%d_c60.E', model_size));
b = B(:,2);
c = C(6,:);
c = c';

r = 3;
color = ['r';'b';'g'];
n = length(A)
rtol = 1e-6;            % RBiCG Tol
itol = 1e-6;            % IRKA Tol
rmaxit = 200;           % RBiCG Max Iteration
imaxit = 50;            % IRKA Max Iteration
irka_iter = 1;

sigma = logspace(-5,0.7,r)';        % Any ways sorts ascending.
fprintf(1,'Real sigma = %14.4e \n', real(sigma));
fprintf(1,'\n');
fprintf(1,'Imag sigma = %14.4e \n', imag(sigma));

for count = 1:r
    x0{count} = zeros(n,1);
    x0_tilde{count} = zeros(n,1);    
end

tic; 
if (recycleFlag == 1)        
    U{1} =[];
    U_tilde{1} = [];
end        
while (irka_iter <= imaxit)

    if (recycleFlag == 1)        
        [V W U U_tilde x0 x0_tilde] = ...
            computeVWRecycle(r, A, E, sigma, droptol, b, c, rtol, ...
            rmaxit, s, k, U, U_tilde, x0, x0_tilde, irka_iter);
    else
        [V W x0 x0_tilde] = ...
            computeVWPlain(r, A, E, sigma, droptol, b, c, rtol, rmaxit, ...
            x0, x0_tilde);
    end
    
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
toc;
disp('All done');               % display the answer.


%% subfunction computeWVRecycle
function [VOrth WOrth U U_tilde x0 x0_tilde] = ...
    computeVWRecycle(r, A, E, sigma, droptol, b, c, tol, maxit, ...
    s, k, U, U_tilde, x0, x0_tilde, irka_iter)

V = [];
W = [];
for i=1:r
    Mat_A = sigma(i)*E - A;
    [Lp,Up] = luinc(Mat_A, droptol);
    bPrecond = Lp\b;
    cPrecond = Up'\c;
    
    if (i == 1)
        if (mod(irka_iter-1, 5) == 0)    % Build space for every 5th system
            [x, x_tilde, ~, ~, resvec, resvec_tilde, U1{i}, U1_tilde{i}] ...
                = rbicg(Mat_A, bPrecond, cPrecond, tol, maxit, x0{i}, ...
                x0_tilde{i}, U{i}, U_tilde{i}, 1, s, k, Lp, Up);
            if isempty(U1{i})
                % do nothing
            else
                U{i} = U1{i};
                U_tilde{i} = U1_tilde{i};
            end
        else
            [x, x_tilde, ~, ~, resvec, resvec_tilde] ...
                = rbicg(Mat_A, bPrecond, cPrecond, tol, maxit, x0{i}, ...
                x0_tilde{i}, U{i}, U_tilde{i}, 0, s, k, Lp, Up);
        end
    else
        [x, x_tilde, ~, ~, resvec, resvec_tilde] ...
            = rbicg(Mat_A, bPrecond, cPrecond, tol, maxit, x0{i}, ...
            x0_tilde{i}, [], [], 0, 0, 0, Lp, Up);
    end    

    length(resvec)
    length(resvec_tilde)
    
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


%% subfunction computeWVPlain
function [VOrth WOrth x0 x0_tilde] = ...
    computeVWPlain(r, A, E, sigma, droptol, b, c, tol, maxit, x0, x0_tilde)

V = [];
W = [];
for i=1:r
    Mat_A = sigma(i)*E - A;
    [Lp,Up] = luinc(Mat_A, droptol);
    bPrecond = Lp\b;
    cPrecond = Up'\c;
    
    [x, x_tilde, ~, ~, resvec, resvec_tilde] = rbicg(Mat_A, bPrecond, ...
        cPrecond, tol, maxit, x0{i}, x0_tilde{i}, [], [], 0, 0, 0, Lp, Up);    

    length(resvec)
    length(resvec_tilde)    

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
