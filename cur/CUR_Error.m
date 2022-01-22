function E = CUR_Error(A,C,R)
%   Error of CUR for given skeleton
%   Independent of the conditioning of U
%   A = [m*n *sparse matrix]
%   C = [m*k *sparse]
%   R = [l*n *sparse]
%   U = argmin_{U} norm(A - C*U*R)
%   E = A - C*U*R
%%
    U = CUR_U_eval(A,C,R);
    E = A - C * (U * R);
end