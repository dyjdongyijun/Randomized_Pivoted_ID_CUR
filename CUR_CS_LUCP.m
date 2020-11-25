function [I,J,U] = CUR_CS_LUCP(A,k,streaming)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%   Outputs:
%       ik, jk = [1*k row / column indices]
%       I = [1*m row idx]; J = [1*n col idx]
%       C = A(:,j(1:k)); R = A(i(1:k),:)
    os = 10;
    [m,n] = size(A);

    if nargin > 2 && streaming
        S_col = embed(m, k+os, 'sparse3');
        S_row = embed(n, k+os, 'sparse3');
        X = S_col(A)'; %(n,k+os)
        Y = S_row(A')'; %(m,k+os)
%         [I,~,~,~,~] = lucp(Y,k);
%         [J,~,~,~,~] = lucp(X,k);
        [I,~] = lucp_lapack(full(Y),k);
        [J,~] = lucp_lapack(full(X),k);
        
    else
        S = embed(n, k+os, 'sparse3');
        Y = S(A')'; %(m,k+os)
%         [I,~,~,~,~] = lucp(Y,k);
%         [J,~,~,~,~] = lucp(A(I(1:k),:)',k);
        [I,~,] = lucp_lapack(full(Y),k);
        [J,~,] = lucp_lapack(full(A(I(1:k),:).'),k);
    end
    
    I = I(1:k); J = J(1:k);
    if nargout==3
        C = A(:,J(1:k)); R = A(I(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end

%% 