function [D] = DifferentiationMatrix(n, xnodes)
    % Input: n - number of nodes 
    % Output: D - differentiation matrix (size n by n)
    %         x -  nodes (size n, column vector)
    

    % Initialize differentiation matrix
    D = zeros(n, n);

    for k=1:n
        xk = xnodes(k);
        
        for i=1:n
            xi = xnodes(i);

            denom = 1; nume = 1;
            for j=1:n
                if j~=k
                    denom = denom*(xk-xnodes(j)) ;
                end
                if j~=k && j~=i
                    nume = nume*(xi-xnodes(j)) ;
                end
            end
            if k~=i
                D(i, k) = nume/denom;
            end
        end
    end

    for i = 1:n
         D(i,i) = -sum(D(i,:));
    end
end




% function [D] = DifferentiationMatrix(n, xnodes)
%     % Input: n - number of nodes 
%     % Output: D - differentiation matrix (size n by n)
%     %         x -  nodes (size n, column vector)
% 
% 
%     % Initialize differentiation matrix
%     D = zeros(n, n);
% 
%     for i=1:n
%         %x = xnodes(i);
%         xnodes(i);
%         % Loop through each Lagrange basis polynomial
%         for k = 1:n
%             xk = xnodes(k); % Current node lagrange basis
% 
%             if k==i
%                 sum = 0;
%                 for j = 1:n
%                     if j ~= k  % Skip the term where j == k
%                         sum = sum + 1/(xk - xnodes(j));
%                     end
%                 end            
%                 D(i,k) = sum;
%             else        
% 
%             numerator = 1;  % Initialize numerator to 1
%             denominator = 1;  % Initialize denominator to 1
% 
%             % Numerator loop (for xnodes(i) - [xnodes(1:k-1)
%             % xnodes(k+1:end)]) skipping one term for addition
%             for j = 1:n
%                 if j ~= k && j~=i % Skip the term where j == k and j==p
%                     numerator = numerator * (xnodes(i) - xnodes(j));
%                 end
%             end
% 
%             % Denominator loop (for xk - [xnodes(1:k-1) xnodes(k+1:end)])
%             for j = 1:n
%                 if j ~= k  % Skip the term where j == k
%                     denominator = denominator * (xk - xnodes(j));
%                 end
%             end
% 
%             % Now calculate the result
%             D(i,k) = (numerator / denominator)  ; 
%             end
%         end
%     end
% 
% end
% 
% 


