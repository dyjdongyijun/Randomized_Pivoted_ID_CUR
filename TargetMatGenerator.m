function targetmat = TargetMatGenerator(name, varargin)
%   {name, *para} = {{'gauss', m, n, k, f},
%                    {'hilbert', m, n},
%                    (disabled){'lapkernel', ns, nt, kernel, helmholtz, plot},
%                    {'gspd',m,k,f}
%                    {'snn',m,n,k,r,a,b,s}: sparse non-negative mat:
%                       A = (m*n), r = rank, k = spectrum threshold, a,b =
%                       multiplier before/after threshold
%                       s = sparsity
%                    {'laplacian',n,m}: graph laplacian of a graph:
%                       G=(V,E), n=|V|, m=|E|, nnz(L_G)=2m+n 
%                   }
%   targetmat = [struct] fieldnames = {
%               'A' = m*n matrix,
%               'description' = matrix description + size
%               'r' = numerical rank: sigma(r+1)<eps*sigma(1)<=sigma(r)
%               'sigma' = [min(m,n)*1] singular values
%               'U','V' = left/right singular vectors
    if nargin < 1, name = 'gauss'; end
%%
    if size(name,1) > 1
        A = name; [m,n] = size(A); 
        if ~isempty(varargin) && strcmpi(varargin{1}, 'svd')
            [U,S,V] = svd(full(A),'econ');
            sigma = diag(S);
        else
            sigma = svd(full(A),'econ');
            U = 0; V = 0;
        end
        k = length(sigma(sigma>eps*sigma(1))); 
        description = sprintf('[Input matrix $%d \\times %d$]',m,n);
%%
    elseif strcmpi(name,'gauss')
        if length(varargin) < 4
            m = 1000; n = 1000; k = 500; f = @(x)exp(-x./2);
        else
            m = varargin{1}; n = varargin{2}; k = varargin{3}; f = varargin{4};
        end
        [A,k,sigma,U,V] = GaussGenerator(m,n,k,f);
        if length(f)==1
            description = sprintf('[Gauss $%d \\times %d$: $%s$]',m,n,func2str(f));
        else
            description = sprintf('[Gauss $%d \\times %d$]',m,n);
        end
%%
    elseif strcmpi(name,'hilbert')
        if length(varargin) < 1
            m = 1000; n = m;
        elseif length(varargin) < 2
            m = varargin{1}; n = m;
        else
            m = varargin{1}; n = varargin{2};
        end
        A = hilb(max(m,n)); A = A(1:m,1:n);
        description = sprintf('[Hilbert $%d \\times %d$]',m,n);
        [U,S,V] = svd(A,'econ');
        sigma = diag(S);
        k = length(sigma(sigma>eps*sigma(1)));
%%
%     elseif strcmpi(name,'lapkernel')
%         [A,k,sigma] = LaplaceKernelMat(varargin{:});
%         [m,n] = size(A);
%         description = sprintf('[Laplace kernel $%d \\times %d$]',m,n);
%         [U,~,V] = svd(A,'econ');
%%
    elseif strcmpi(name,'gspd')
        if length(varargin) < 3
            m = 1000; k = 200; f = @(x)exp(-x./2);
        else
            m = varargin{1}; k = varargin{2}; f = varargin{3};
        end
        [A,k,sigma,U] = GaussSPDGenerator(m,k,f); 
        V = U;
        if length(f)==1
            description = sprintf('[Gauss SPD $%d \\times %d$: $%s$]',m,m,func2str(f));
        else
            description = sprintf('[Gauss SPD $%d \\times %d$]',m,m);
        end
%%
    elseif strcmpi(name,'snn')
        if length(varargin) < 7
            s = 2e-3;
        else
            s = varargin{7};
        end
        if length(varargin) < 2
            m = 300000; n = 300;
        else
            m = varargin{1}; n = varargin{2};
            if  length(varargin) < 6
                k = min([10,m,n]); r = min([m,n]); a = 2; b = 1;
            else
                k = varargin{3}; r = varargin{4}; 
                a = varargin{5}; b = varargin{6};
            end
        end
        X = sprand(m,r,s); Y = sprand(r,n,s); 
        d = 1./(1:r);
        d(1:k) = d(1:k)*a; d(k+1:r) = d(k+1:r)*b;
        A = (X.*d)*Y;
        description = sprintf('[Sparse non-negative $%d \\times %d$, $nnz=%d$]',m,n,nnz(A));
        [U,S,V,flag] = svds(A, min( r, max(floor(min(m,n)/10),5000) ));
        if flag==1 && (m+n)<2e4 
            [U,S,V] = svd(full(A),0);
        end
        sigma = diag(S);
        k = length(sigma(sigma>eps*sigma(1)));
%%
    elseif strcmpi(name,'laplacian')
        if length(varargin) < 2
            n = 1000; m = 5000;
        else
            n = varargin{1}; m = varargin{2};
        end
        A = rand_laplacian(n,m);
        description = sprintf('[Random Laplacian: $n=%d$, $m=%d$]',n,m);
        [U,S,V,flag] = svds(A, max(floor(n/10), 5000));
        if flag==1 && (m+n)<2e4 
            warning('svds does not converge')
            [U,S,V] = svd(full(A),0);
        end
        sigma = diag(S);
        k = length(sigma(sigma>eps*sigma(1)));
    end
%%
    targetmat = struct('A',{A},...
        'description',{description},...
        'r',{k},...
        'sigma',{sigma});
    targetmat.('U') = U;
    targetmat.('V') = V;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,k,s,U,V] = GaussGenerator(m,n,k,f)
%   generate m*n random Gaussion matrix A with rank k via SVD 
%   f = decay function of singular values: d(1:k) = sort(f(1:k),'descend')
%   or  1*length array with length >= k
    [U,~] = qr(randn(m,k),0);
    [V,~] = qr(randn(n,k),0);
    d = sort(abs(f(1:k)),'descend');
    A = U*(d.*V)';
    k = length(d(d>eps*d(1)));
    s = zeros(1,min(m,n));
    s(1:length(d)) = d;
end

function [A,k,s,U] = GaussSPDGenerator(m,k,f)
%   generate m*m random Gaussion SPD matrix A with rank k via spectral decomp 
%   f = decay function of singular values: d(1:k) = sort(f(1:k),'descend')
%   or  1*length array with length >= k
    [U,~] = qr(randn(m,k),0);
    d = sort(abs(f(1:k)),'descend');
    A = U*(d.*U)';
    k = length(d(d>eps*d(1)));
    s = zeros(1,m);
    s(1:length(d)) = d;
end
