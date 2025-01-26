function [lambda V V_Res W W_Res] = getGenEigenvecs(k, A, B)
% We need both the left and the right eigenvectors, and thus we use "qz"
% function. By default it does a complex decomposition with the "AA" and
% "BB", while we need all real.
[AA, BB, Q, Z, V, W] = qz(A,B,'real');

% Sort the eigenvalues.
% - abs ensures the magnitude is picked for both negative and complex val.
% - QZ is fine in  maintaining the order. See the Math Res II. 
lambda = ordeig(AA,BB);
[lambda_sorted, lambda_idx] = sort(abs(lambda));
V = V(:, lambda_idx);
W = W(:, lambda_idx);
lambda = lambda(lambda_idx);

% Pick the smallest (absolute) k. Also, making sure that complex conjugate
% pairs are carried forward. Note that as we are working with real matrics
% so both the eigenvectors and eigenvalues occur in complex conjugate
% pairs. Also, when comparing eigenvectors one has to worry of the
% normalization, but not when working with eigenvalues.
ev1 = lambda(k);
ev2 = lambda(k+1);
if abs(real(ev1) - real(ev2)) <= 1e-06 && abs(imag(ev1) - imag(ev2)) <= 1e-06
    endIndex = k+1;
else
    endIndex = k;
end
V = V(:, 1:endIndex);
W = W(:, 1:endIndex);
lambda = diag(lambda(1:endIndex));

V_Res = [];
W_Res = [];

% % Prepare the V and W such that the residual: A*V_Res - V_Res*lambda could
% % be computed.
% count = 1;
% while (count <= endIndex)
%     if (count <= endIndex-1) && (AA(count+1, count) ~= 0)
%         V_Res(:,count) = V(:,count) + i*V(:,count+1);
%         W_Res(:,count) = W(:,count) + i*W(:,count+1);
%         V_Res(:,count+1) = V(:,count) - i*V(:,count+1);
%         W_Res(:,count+1) = W(:,count) - i*W(:,count+1);
%         count = count + 1;
%     else
%         V_Res(:,count) = V(:,count);
%         W_Res(:,count) = W(:,count);
%     end
%     count = count + 1;
% end

