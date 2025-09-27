function w = quad_weights_lobatto(tau)
% quad_weights_lobatto: Compute weights for Legendre-Gauss-Lobatto nodes
% Input: tau = [tau_0, ..., tau_n] — LGL nodes (including ±1)
% Output: w = [w_0, ..., w_n] — Lobatto quadrature weights

    n = length(tau);            % total n nodes → L_{n-1}(τ)
    w = zeros(size(tau));
    
    w(1)   = 2 / (n*(n-1));     % at τ_0 = -1
    w(end) = 2 / (n*(n-1));     % at τ_n = 1

    % Interior weights using L_{n-1}(tau)
    for j = 2:n-1
        Lj_val = legendreP(n-1, tau(j));
        w(j) = 2 / (n*(n-1) * Lj_val^2);
    end
end
