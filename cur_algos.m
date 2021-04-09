function [i,j] = cur_algos(algo, A, l)
%   {'SRCUR',...
%    'LUPP', 'LUPP2pass', 'CSLUPP', 'LUPPstream',...
%    'CPQR', 'CPQR2pass', 'CSCPQR', 'CPQRstream',...
%    'RSVDDEIM',...
%    'RSVDLS'}
    switch algo
        % srcur
        case 'SRCUR'
            [i,j] = srcur(A, l);
        % lupp
        case 'LUPP'
            sketch = 'gauss';
            Smpr = embed(m, l, sketch);
            X = Smpr(A); % (l,n)
            [~,~,j] = lu(X','vector');
            C = A(:,j(1:l)); %(m,l)
            [~,~,i] = lu(C,'vector');
        case 'LUPP2pass'
            sketch = 'gauss';
            Smpr = embed(n, l, sketch);
            Yt = Smpr(A'); % (l,m)
            X = Yt*A; %(l,n)
            [~,~,i] = lu(X','vector');
            R = A(i(1:l),:); %(l,n)
            [~,~,j] = lu(R','vector');
        case 'CSLUPP'
            sketch = 'sparse8';
            Smpr = embed(m, l, sketch);
            X = Smpr(A); % (l,n)
            [~,~,j] = lu(X','vector');
            C = A(:,j(1:l)); %(m,l)
            [~,~,i] = lu(C,'vector');
        case 'LUPPstream'
            sketch = 'gauss';
            Smpr_m = embed(m, l, sketch);
            Smpr_n = embed(n, l, sketch);
            X = Smpr_m(A); % (l,n)
            Yt = Smpr_n(A'); % (l,m)
            [~,~,j] = lu(X','vector');
            [~,~,i] = lu(Yt','vector');
        % cpqr
        case 'CPQR'
            sketch = 'gauss';
            Smpr = embed(m, l, sketch);
            X = Smpr(A); % (l,n)
            [~,~,j] = qr(full(X),0);
            C = A(:,j(1:l)); %(m,l)
            [~,~,i] = qr(full(C'),0);
        case 'CPQR2pass'
            sketch = 'gauss';
            Smpr = embed(n, l, sketch);
            Yt = Smpr(A'); % (l,m)
            X = Yt*A; %(l,n)
            [~,~,i] = qr(full(X),0);
            R = A(i(1:l),:); %(l,n)
            [~,~,j] = qr(full(R),0);
        case 'CSCPQR'
            sketch = 'sparse8';
            Smpr = embed(m, l, sketch);
            X = Smpr(A); % (l,n)
            [~,~,j] = qr(full(X),0);
            C = A(:,j(1:l)); %(m,l)
            [~,~,i] = qr(full(C'),0);
        case 'CPQRstream'
            sketch = 'gauss';
            Smpr_m = embed(m, l, sketch);
            Smpr_n = embed(n, l, sketch);
            X = Smpr_m(A); % (l,n)
            Yt = Smpr_n(A'); % (l,m)
            [~,~,j] = qr(full(X),0);
            [~,~,i] = qr(full(Yt),0);
        % deim
        case 'RSVDDEIM'
            [V,~,W] = RSVD(A,l);
            [~,~,i] = lu(V,'vector');
            [~,~,j] = lu(W,'vector');
        % leverage score sampling
        case 'RSVDLS'
            [V,~,W] = RSVD(A,l);
            i = LeverageScoreSampling(V); 
            j = LeverageScoreSampling(W); 
        otherwise % using LUPP
            sketch = 'gauss';
            Smpr = embed(m, l, sketch);
            X = Smpr(A); % (l,n)
            [~,~,j] = lu(X','vector');
            C = A(:,j(1:l)); %(m,l)
            [~,~,i] = lu(C,'vector');
    end    
end

function p = LeverageScoreSampling(V)
%   V = [m*k matrix] top k left / right singular vectors
%   p = [1*k] row skeleton indices selected by LS
    ls = sum(V.^2,2)./size(V,2);
%     p = sampleProb((1:size(V,1))', ls, size(V,2));
    [~,p] = maxk(ls, size(V,2));
end