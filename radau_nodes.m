function nodes = radau_nodes(n)
% Rn(τ) = L̃n(τ) − L̃n−1(τ),     where τ ∈ [−1, 1]
% radau_nodes  Roots of flipped Legendre‑Radau polynomial R_n(x) 
% Returns (n) flipped Legendre‑Radau nodes in [-1, 1]
syms x;
roots = vpasolve(legendreP(n+1,x) - legendreP(n,x) == 0);
nodes = sort(roots);
end

% function nodes = radau_nodes(n)
% % radau_nodes  Roots of flipped Legendre‑Radau polynomial R_n(x)
% 
%     Ln   = legendreP_coeffs(n);
%     Ln_1 = legendreP_coeffs(n-1);
%     Rn   = Ln - [0 Ln_1];          % degree n
%     rts  = roots(Rn);              % complex in general
%     tol  = 1e-12;
%     rts  = rts(abs(imag(rts))<tol);% keep real
%     nodes = sort(real(rts));
% end

% function x = radau_nodes(n)
% % radau_nodes  Left-Radau nodes: n nodes in [-1,1], includes -1, excludes +1
% % Uses roots of Jacobi polynomial P_{n-1}^{(0,1)}(x)
% 
%     if n < 2
%         error('n must be >= 2 for Radau nodes');
%     end
% 
%     m = n - 1;           % Number of interior nodes (excluding -1)
% 
%     % Jacobi matrix for P_{m}^{(0,1)}(x)
%     i = (1:m-1)';
%     a = zeros(m,1);  % Diagonal entries of Jacobi matrix
%     b = sqrt(i .* (i + 1) ./ ((2*i + 1) .* (2*i + 3)));  % Off-diagonal
% 
%     J = diag(a) + diag(b,1) + diag(b,-1);
% 
%     % Compute interior nodes (roots in (-1, 1))
%     interior_nodes = sort(eig(J));
% 
%     % Left-Radau: include -1
%     x = sort([-1; interior_nodes]);
% end
