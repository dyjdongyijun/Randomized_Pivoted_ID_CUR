function [I, J] = lucp_lapack(A,k)
    if ~exist('k','var')
        k = min(size(A));
    end
    
    if size(A,1) >= size(A,2)
        LDA = size(A,1);
        N = size(A,2);
        I = zeros(N,1);
        J = zeros(N,1);
        flag = -1;
        [~,~,~,I_lpk,J_lpk,flag] = lapack('dgetc2',N,A,LDA,I,J,flag);
        I = lapack2perm(I_lpk, k);
        J = lapack2perm(J_lpk, k);

    else
        LDA = size(A,2);
        N = size(A,1);
        I = zeros(N,1);
        J = zeros(N,1);
        flag = -1;
        [~,~,~,J_lpk,I_lpk,flag] = lapack('dgetc2',N,A.',LDA,J,I,flag);
        I = lapack2perm(I_lpk, k);
        J = lapack2perm(J_lpk, k);
        
    end
    
    if flag~=0
        warning('U(%d, %d) is likely to produce overflow if we try to solve for x in Ax = b. So U is perturbed to avoid the overflow.',...
                flag, flag)
    end
end

function I = lapack2perm(Ilpk, k)
    I = 1:length(Ilpk);
    for i = 1:k
        aux = I(Ilpk(i));
        I(Ilpk(i)) = I(i);
        I(i) = aux;
    end
end