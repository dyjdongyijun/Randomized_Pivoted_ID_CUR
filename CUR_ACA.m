function [I,J,U] = CUR_ACA(A,k)
%   A = [m*n matrix]
%   k = target rank
%   I= [1*k subset of [m]] row skeletons
%   J= [1*k subset of [n]] column skeletons
%   U = [k*k matrix]

    os = 10; 
    maxiter = 5;    
%     p = floor((min(k*log(m)/eps^2, m)));
%     q = floor(min(k*log(n)/eps^2, n));
    [m,n] = size(A);
    l = k + os;
%%  row skeleton
    S_row = embed(n,l,'gauss');
    Q = S_row(A')'; % (m,l)
    [Iout,~,~,~,kout] = adaptive_cross_approx(Q,k);
    i = 0;
    while kout < k && i < maxiter
        [Iout,~,~,~,kout] = adaptive_cross_approx(Q,k);
        i = i+1;
    end
    I = Iout(1:k);
    
%%  column skeleton
    S_col = embed(m,l,'gauss');
    P = S_col(A)'; % (n,l)
    [Jout,~,~,~,kout] = adaptive_cross_approx(P,k);
    i = 0;
    while kout < k && i < maxiter
        [Jout,~,~,~,kout] = adaptive_cross_approx(P,k);
        i = i+1;
    end
    J = Jout(1:k);
%%
    I = I(1:k); J = J(1:k);
    if nargout==3
        C = A(:,J(1:k)); R = A(I(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end