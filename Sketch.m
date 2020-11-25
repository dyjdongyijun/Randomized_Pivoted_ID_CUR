function Y = Sketch(A,k,randmat)
%   inputs: A = (m*n) matrix
%           k = sampling size, i.e., number of columns of the output
%           randmat = {'gauss','fsrft','srft','gsrft','sparse(n)'}
%   output: Y = (m*k) (structured) sketch of the column space of A
%%
    if nargin < 3, randmat = 'gauss'; end
    if strcmpi(randmat,'fsrft')
       Y = SRFT(A,k);
    elseif strcmpi(randmat,'srft')
       Y = SRFT(A,k,'srft');
    elseif strcmpi(randmat,'gsrft')
       Y = GSRFT(A,k);
    elseif length(randmat)>=length('sparse') && strcmpi(randmat(1:length('sparse')),'sparse')
       if length(randmat)==length('sparse')
           s = 3;
       else
           s = str2double(randmat(length('sparse')+1:length(randmat)));
       end
       Y = RandSparse(A,k,s);
    else % if strcmpi(randmat,'gauss')
       G = randn(size(A,2),k);
       Y = A*G;
    end
end

%%
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

%%
function Y = GSRFT(A,l)
%   Compute sampling of A by random chain of Givens rotations: 
%   B = A*D2*T2*D1*T1*D*F*S
    Y = A;
    for i = 1:2
        d = exp((1i*2*pi).*rand(1,size(A,2)));
        Y = RandGivensChain(d.*Y);
    end
    Y = SRFT(Y,l,'full');
end

%% %%%%%%%%%%% Theta = Chain of random Givens rotation %%%%%%%%%%%%%
function B = RandGivensChain(A)
%   B = A \Theta, where \Theta = \Pi*G(1,2;\theta_1)*...*G(n-1,n;\theta_{n-1})
    pcol = randperm(size(A,2));
%   Apply B=A*\Pi
    B = A(:,pcol);
%   Apply chain of random Givens rotations
    thetas = rand(1,size(A,2)-1).*(2*pi);
    cs = cos(thetas); ss = sin(thetas);
%   Can we get rid of the loop? (No since the rotations are not commutative?)
    for i = 1:size(A,2)-1
       B(:,i:i+1) = B(:,i:i+1)*[cs(i), -ss(i); ss(i), cs(i)]; 
    end
end

%%
function Y = RandSparse(A,l,s)
%   Compute sampling of A by random sparse matrix: B = A*S
%   A = (m*n), rank(A) = r < l
%   l = sampling size, if l<s, l=s; end
%   Y = (m*l) sample matrix
%   s = number of nonzero entries in each row of S
%       total number of nonzeros = s*n
    i = repelem(1:size(A,2),s);
    j = zeros(size(i));
    for row = 1:size(A,2)
        j((row-1)*s+1 : row*s) = randperm(max(l, s), s);
    end
%     j = cell2mat(arrayfun(@(aux) randperm(max(l,s),s), 1:size(A,2),...
%                         'UniformOutput', false));
    S = sparse(i,j,(-1).^randi(2,size(i)),size(A,2),max(l,s));
    Y = A*S(:,1:l);
end
