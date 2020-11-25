function [I,J,t,U,V] = lucp(A,k,compute_uv)
%   Inputs: A = [m*n matrix], 
%           k = target rank, 
%   Outputs: I, J = [1*m, 1*n index vector] row/col indices selected in top-t entries,
%            U = [m*t] unit lower triangular, 
%            V = [n*t] upper triangular,
%            A(I,J) \approx U*V'

    eps = 1e-10;
    
    [m,n] = size(A);
    k = min([n,m,k]);
    I = 1:m; J = 1:n;
    Iact = I; Jact = J;
    
    for t = 1:k
%       global search for the max
        [colmaxval, colmaxidx] = max(abs(A(Iact,Jact)));
        [pivot,jaux] = max(colmaxval);
        if pivot < eps
            t = t-1;
            break;
        end
        iaux = colmaxidx(jaux);
        ip = Iact(iaux); jp = Jact(jaux);
%       swap indices
        I(iaux+t-1) = I(t);
        I(t) = ip;
        J(jaux+t-1) = J(t);
        J(t) = jp;
%       update active range
        Iact = I(t+1:end);
        Jact = J(t+1:end);
%       update pivoted subcolumn in-place
        u = A(Iact,jp)./A(ip,jp);
        A(Iact,jp) = u;
%       update submatrix in-place
        v = A(ip,Jact);
        A(Iact,Jact) = A(Iact,Jact) - u*v;
    end
    
    if nargin >= 3 && compute_uv
        U = tril(A(I,J(1:t)),-1)+speye(m,t);
        V = triu(A(I(1:t),J))';
    else
        U = 0; V = 0;
    end
end