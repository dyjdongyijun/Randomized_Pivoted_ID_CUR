function LG = rand_laplacian(n,m)
%   Inputs
%       n: [int] number of vertices2E
%       m: [int] number of edges
%   Output: G with weights uniform [0,1]
%       LG: [n*n sparse] Laplacian of G=(V,E), |V|=n, |E|=m, nnz(LG)=2*m+n
    adj = spalloc(n, n, m);
    smp = 0;
    while smp < m
        smp = min(max(3*m, smp+m), floor(n*n/2));
        idx = randperm(n*n, smp);
        i = ceil(idx/n);
        j = mod(idx-1,n)+1;
        idx(i>=j) = []; % remove diagonal and upper triangluar indices
        smp = length(idx);
    end
    adj(idx(1:m)) = rand(1,m);
    adj = adj + adj';
    dia = diag(sum(adj));
    LG = dia - adj;
end