function [i,j,U] = CUR_ID_streaming(A,k,randmat,l,s)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%       randmat = {'gauss','srft','sparse(n)'}
%       l = row sample size: m -> l
%       s = col sample size: n -> s
%   Outputs:
%       ik, jk = [1*k row / column indices]
%       i = [1*m row idx]; j = [1*n col idx]
%       C = A(:,j(1:k)); R = A(i(1:k),:) 
%       U = [k*k matrix]
    if nargin < 3, randmat='gauss'; end
    if nargin < 4, l = min(min(size(A)),8*k); s = l; end
    if nargin < 5, s = l; end
    
    [m,n] = size(A);
    S_col = embed(m, l, randmat);
    S_row = embed(n, s, randmat);
    X = S_col(A)'; % (n,l)
    Y = S_row(A')'; % (m,s)
    [~,~,j] = qr(full(X'),0);
    [~,~,i] = qr(full(Y'),0);
    if nargout == 3
        Z = ( S_col(Y) + S_row(X)' )./2; % (l,s)
                
        Chat = X(j(1:k),:)'; %(l,k)
        Rhat = Y(i(1:k),:); %(k,s);
        
        K = Chat\Z;
        U = (Rhat'\K')';
    end
        
end