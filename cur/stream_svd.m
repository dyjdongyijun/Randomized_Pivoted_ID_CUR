function [U,d,V] = stream_svd(A, k, randmat, l, s)
% Input
%   A = [m*n matrix]
%   k << min(m,n) sample rank
%   l = sample size for left / right sketching
%   s = sample size for core sketching 
%       (if s <= 0, no core sketch)
%   randmat = {'srft', 'sparse(n)'}
% Output
%   U = [m*k]
%   s = [k*k]
%   V = [n*k]

    if nargin<3, randmat = 'sparse3'; end
    if nargin<4, l = 2*k; end
    if nargin<5, s = 2*l; end
    
    if s > 0
        s = ceil(max(s, 1.1*l));
        [U,d,V] = stream_svd_3_sketch(A, k, randmat, l, s);
    else
        [U,d,V] = stream_svd_2_sketch(A, k, randmat, l);
    end
end


function [U,d,V] = stream_svd_3_sketch(A, k, randmat, l, s)
    [m,n] = size(A);

    S_row = embed(m, l, randmat);
    S_col = embed(n, l, randmat);
    S_core_left = embed(m, s, randmat);
    S_core_right = embed(n, s, randmat);
    
    X = S_row(A)';  % (n,l)
    Y = S_col(A')'; % (m,l)
    Z = S_core_right(S_core_left(A)')'; % (s,s)
    
    [P,~] = qr(X,0); % (n,l)
    [Q,~] = qr(Y,0); % (m,l)
    
    coef_left = S_core_left(Q); % (s,l)
    coef_right = S_core_right(P)'; % (l,s)
    
    solve_left = coef_left\Z; % (l,s)
    C = (coef_right'\solve_left')'; % (l,l)
    
    [Ul, Dl, Vl] = svd(full(C));
    U = Q*Ul(:,1:k);
    V = P*Vl(:,1:k);
    dl = diag(Dl); d = dl(1:k);
end


function [U,d,V] = stream_svd_2_sketch(A, k, randmat, l)
    [m,n] = size(A);

    S_row = embed(m, l, randmat);
    S_col = embed(n, l, randmat);
    
    X = S_row(A)';  % (n,l)
    Y = S_col(A')'; % (m,l)
    
    [P,~] = qr(X,0); % (n,l)
    [Q,~] = qr(Y,0); % (m,l)
    
    coef_right = S_col(P)'; % (l,l)
    C = (coef_right'\(Q'*Y)')'; % (l,l)
    
    [Ul, Dl, Vl] = svd(full(C));
    U = Q*Ul(:,1:k);
    V = P*Vl(:,1:k);
    dl = diag(Dl); d = dl(1:k);
end