function w = quad_weights_radau(tau)
% quad_weights_radau: Compute Radau quadrature weights using flipped Radau nodes
% Input:  tau = [tau_1, ..., tau_n] — roots of flipped Legendre-Radau poly
% Output: w   = [w_1, ..., w_n] — Radau quadrature weights
%
% Note: tau(1) should be -1 (Radau endpoint); rest are interior nodes in [-1,1].

    n = length(tau);
    %rev_tau = -sort(-tau);

    syms t

    % Define Legendre polynomials symbolically
    L_n   = legendreP(n, t);
    
    % Convert to function handles
    L_n_fun = matlabFunction(L_n, 'Vars', t);

    w_tilde = zeros(size(tau));
    
    % % First weight: analytical formula for τ = -1
    w_tilde(1) = 2 / n^2;
    
    % Compute remaining weights
    for j = 2:n
        %tau_j = rev_tau(j);
        tau_j = tau(j);
        Lj_val = L_n_fun(tau_j);
        w_tilde(j) = (1 - tau_j) / (n^2 * Lj_val^2);
    end
    
    % Flip operator (for flipped Radau): multiply by -1 and sort
    w = abs(sort(-w_tilde));
end
