function nodes = lobatto_nodes(n)
% Ln(τ) = (1 - τ^2) * d/dτ [Ln−1(τ)],     where τ ∈ [−1, 1]
% lobatto_nodes  Legendre-Lobatto nodes (n points)
% Returns (n) LGL nodes in [-1, 1]
syms x;
roots = double(vpasolve(diff(legendreP(n-1,x)) == 0));
nodes = sort([-1; roots; 1]);
end


% function nodes = lobatto_nodes(n)
% % lobatto_nodes  Legendre-Lobatto nodes (n+1 points)
% % Returns (n+1) LGL nodes in [-1, 1]
% 
%     if n < 1, error('n must be ≥ 1'); end
% 
%     % Derivative of Legendre polynomial of degree n
%     Ln = legendreP_coeffs(n);       % Get coefficients of P_n(τ)
%     dLn = polyder(Ln);              % Derivative: P_n'(τ)
% 
%     % Roots of the derivative (interior Lobatto nodes)
%     rts = roots(dLn);
%     tol = 1e-12;
%     %rts = rts(abs(rts) <= 1 + 1e-12); % discard numerical drift outside [-1,1]
%     rts = rts(abs(imag(rts)) < tol);  % Drop imaginary roots
%     rts = real(rts);
% 
%     % Add endpoints
%     nodes = sort([-1; rts; 1]);
% end

% function x = lobatto_nodes(n)
%     % lobatto_nodes returns n+1 Legendre-Gauss-Lobatto nodes in [-1,1]
%     % Uses the stable Golub-Welsch method for Legendre polynomial roots
% 
%     if n < 1
%         error('n must be >= 1');
%     end
% 
%     % Compute interior Legendre points (roots of P'_{n})
%     % Use the Jacobi matrix for stable computation
%     i = 1:n-1;
%     a = zeros(n-1,1);
%     b = sqrt(i.^2 ./ (4*i.^2 - 1));
%     J = diag(b,1) + diag(b,-1);
%     interior_nodes = sort(eig(J));  % These are in (-1,1)
% 
%     % Add endpoints -1 and 1
%     x = [-1; interior_nodes; 1];
% end





