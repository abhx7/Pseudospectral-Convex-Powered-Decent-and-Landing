function D = diff_mat_barycentric(x)
% diff_mat_barycentric  Return the (N×N) first‑derivative matrix
%   using the barycentric‑form formula (Fornberg, 1996).
    N = numel(x);
    c = ones(N,1);
    for j = 1:N
        c(j) = 1/prod(x(j) - x([1:j-1 j+1:N]));
    end
    D = zeros(N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i,j) = c(j)/(c(i)*(x(i)-x(j)));
            end
        end
        D(i,i) = -sum(D(i,[1:i-1 i+1:N]));
    end
end
