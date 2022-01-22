function [I,J,U] = CUR_LUCP(A,k,disable_lapack)
%   Inputs:
%       A = [m*n matrix]
%       k = [int] sample size = O(rank(A))
%   Outputs:
%       I = [1*k row idx]; J = [1*k col idx]
%       C = A(:,j(1:k)); R = A(i(1:k),:)
    if exist('disable_lapack','var') && disable_lapack
        [I,J,~,~,~] = lucp(A,k);
    else
        [I,J] = lucp_lapack(full(A),k);
    end

    I = I(1:k); J = J(1:k);
    if nargout==3
        C = A(:,J(1:k)); R = A(I(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end