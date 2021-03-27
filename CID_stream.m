function [J, V, etamaxf] = CID_stream(A, k, pivot)
    m = size(A,1); n = size(A,2);
    S = randn(k,m);
    X = S*A;
    
    if strcmpi(pivot(1:4),'cpqr')
        [~,T,Jall] = qr(X, 'vector');
    else
        [L,~,Jall] = lu(X', 'vector');
        T = L';
    end
    T1 = T(:,1:k);
    T2 = T(:,k+1:end);
    aux = T1\T2;
    VT = zeros(k,n);
    VT(:,Jall) = [eye(k), aux];
    V = VT';
    J = Jall(1:k);
    etamaxf = norm(aux, 'fro');
end