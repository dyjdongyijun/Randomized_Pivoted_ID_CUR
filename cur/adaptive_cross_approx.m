function [I, J, U, V, kout] = adaptive_cross_approx(A,k)
%   Inputs: A = [m*n matrix], 
%           k = target rank, 
%   Outputs: I, J = [1*k index vector] row/col indices selected
%            U = [m*kout] unit lower triangular under permutation, 
%            V = [n*kout] upper triangular under permutation,
%            A \approx U*V'

    eps = 1e-10;

    [m,n] = size(A);
    U = zeros(m,k); V = zeros(n,k);
    I = 1:m; J = 1:n;
    
    j = randi(n); 
    Iact = I; Jact = J;
    
    for t = 1:k
        if ~(j==1)
            J(j+t-1) = J(t);
            J(t) = Jact(j);
        end
        
        u = A(Iact,Jact(j)) - U(Iact,1:t-1)*V(Jact(j),1:t-1)';
        [umax,i] = max(abs(u)); 
        if umax < eps 
            t = t-1;
            break; 
        end
        if ~(i==1)
            I(i+t-1) = I(t);
            I(t) = Iact(i);
        end
        u = u./u(i);
        v = A(Iact(i),Jact)' - V(Jact,1:t-1)*U(Iact(i),1:t-1)';
        
%         dfro = sum(u.^2)*sum(v.^2);
%         fro = fro + dfro + sum((V(Jact,1:t-1)'*v).*(U(Iact,1:t-1)'*u));  
        
        U(Iact,t) = u; 
        V(Jact,t) = v;
        
        Iact = I(t+1:end);
        Jact = J(t+1:end);
        [~,j] = max(abs(V(Jact,t))); 
    end
    
    I = I(1:k); J = J(1:k);
    kout = t;
    U = U(:,1:kout); V = V(:,1:kout);
end

