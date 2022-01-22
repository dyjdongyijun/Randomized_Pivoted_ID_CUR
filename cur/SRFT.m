function Y = SRFT(A,l,opt)
%   Compute SRFT: B = A*(DFS)
%   A = p*n matrix, want the l-point SRFT specified by s of the rows
%   l = sample size (l<n): 1. mod(n,l)==0 -> SRFT 
%                          2. otherwise || opt=='full' -> truncated full FFT
%   B = m*l SRFT output
%   opt = [options] 1. full = using full FFT -> O(n*log(n));
%                   2. (otherwise, only when l divides n) = using SRFT ->
%                   O(n*log(l))
    [p,n] = size(A);
    if nargin < 3, opt='full'; end
    if l > n, l=n; opt='full'; end
    if issparse(A), A = full(A); end
%%  Making sure naug = l*m or Full FFT
    m = ceil(n/l); 
    if l*m > n || strcmpi(opt,'full') % use full FFT
        d = exp(1i.*rand(1,n).*2.*pi);
        s = randperm(n,l);
        Yaug = fft(d.*A,[],2);
        Y = sqrt(size(A,2)/l)*Yaug(:,s);
%         warning('Sample size l does not divide n: full FFT is used!')
        return
    end
%%  SRFT
    d = exp(1i.*rand(1,n).*2.*pi);
    s = randperm(n,l);
%   Step 1
    V = reshape(permute(reshape((d.*A)',m,l,[]),[2,1,3]),l,[]);
    W = fft(V,l,1);
%   Step 2
    X = exp((reshape(0:l-1,[],1)*reshape(mod(0:m*p-1,m),1,[])).*(-2*pi*1i/n)).*W;
%   Step 3
    [j1,j2] = index2coord(s,l);
    Y = permute(reshape(X(j2,:),l,[],p),[2,1,3]); % l*m*p
    F = exp((-2*pi*1i/m).*(reshape(0:m-1,[],1)*reshape(j1-1,1,[])));
    Y = squeeze(sum(F.*Y,1))';
end

function [j1,j2] = index2coord(j,l)
%   output matrix indices as row vectors
    j1 = reshape(floor((j-1)/l)+1,1,[]); 
    j2 = reshape(mod(j-1,l)+1,1,[]);
end