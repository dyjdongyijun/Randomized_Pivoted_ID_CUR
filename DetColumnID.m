function [j,V] = DetColumnID(A,k)
%   Inputs:
%       A = [l*n matrix]
%       k = [int] k <= l < n
%   Outputs: A \approx C*V'
%       j = [n column indices in [n]] C:=A(:,j), first k valid
%       V = [n*k matrix] V(j,:)' = [eye(k), T]
    n = size(A,2);
    [~,R,j] = qr(full(A),0); % R:(l,n), jl:(n)
    V = zeros(n,k);
    V(j,:) = [eye(k), R(1:k,1:k)\R(1:k,k+1:end)]'; % (n,k)
%     j = j(1:k);
end
