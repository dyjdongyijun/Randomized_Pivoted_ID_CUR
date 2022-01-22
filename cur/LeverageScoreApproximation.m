function [lst,r1,r2] = LeverageScoreApproximation(A, eps)
% % % 
% Based on [Drineas2012]: arxiv1109.3843
% Approximate the row leverage scores in o(nd^2) time
% % % 
%   A = [n*d], eps <= 1/2
%   lst = approx leverage scores = [n*1]
    delta = 0.1;
    B = FJLT(A, eps);
    r1 = size(B,1);
%     [Ut,St,Vt] = svd(B,'econ'); 
%     R = Vt./diag(St)';
    [~,R] = qr(B);
    C = R'\A';
    Omg = JLT(C, eps, delta)';
    r2 = size(Omg,2);
    lst = sum(Omg.^2, 2);
end

function Y = JLT(X, eps, delta)
%   X = [d*n] n pts in dim=d
%   Y = [r*n] dim-r subspace embedding of n dim-d vectors in X
%   Y = Pi*X, Pi = [r*d], pi = \pm \sqrt{3/r} w.p. 1/6 * 2, = 0 w.p. 2/3
    [d, n] = size(X);
%     r = (12 * log(n) + 6 * log(1/delta))/eps^2;
    r = log(n/delta)/eps^2;
    r = ceil(min(r, d/(1e3)));
    cnz = randperm(r*d, floor(r*d/3)); % coordinates of nonzeros: d*(i-1)+j
    jnz = mod(cnz-1, d)+1;
    inz = ceil(cnz/d);
    val = (-1).^randi(2,size(cnz));
    Pi = sparse(inz, jnz, val, r, d);
    Y = Pi*X;
end

function B = FJLT(A, eps)
%   A = [n*d]
%   B = [r*d] = S'*Hn*D*A
%       D = [n*n] = diag(\pm 1) w.p. 1/2 * 2
%       Hn = [n*n] Hadamard
%       S = n*r partial permutation
%       Fn = [n*n] fft is used here instead
    [n,d] = size(A);
%     r = 14^2*d*log(40*n*d)/eps^2*log(30^2*d*log(40*n*d)/eps^2);
    r = d*log(n)/eps^2*log(d*log(n)/eps^2);
    r = ceil(min(r, n/(1e3)));
    d = (-1).^randi(2,[n,1]);
    DA = d.*A;
    FDA = abs(fft(full(DA)));
    B = FDA(randperm(n,r),:);
end