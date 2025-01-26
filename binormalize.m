function [C_hat, C_tilde_hat] = binormalize(C, C_tilde, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize columns of C and C_tilde with c1*c1_tilde etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_hat = C_tilde;            % C_tilde with normalized columns.
C_tilde_hat = C;            % C with normalized columns.
for i = 1:k
    normScalar = C(:,i)'*C_tilde(:,i);
    C_hat(:,i) = C_hat(:,i)/normScalar;
    C_tilde_hat(:,i) = C_tilde_hat(:,i)/normScalar;
end