function [J,VT] = rand_id(A,tol,k_max,option,enable_full_cpqr)
% Inputs:
%   A: (m,n)
%   tol: norm(A-A(:,J)*VT) \lesssim tol
%   (k_max): (default=floor(n/log2(n))) over-estimated rank
%   (option): (default='lupp') {'lupp','cpqr'} with sketching
%   (enable_full_cpqr): (default=false) switch to full CPQR if k_max is not 
%       sufficient for the tolerance

% Outputs:
%   J: (k) indices of skeleton columns
%   (V): (n,k) such that A(:,J)*VT approximates A

% Hyperparameters:
%   ETA_L_POW_CPQR (for cpqr), ETA_L_POW_SVD (for LUPP): 
%       Increasing the constants improves the accuracy, while compromises
%       the efficiency.
    
    m = size(A,1);
    n = size(A,2);
    ETA_L_POW_CPQR = 0.5;
    ETA_L_POW_SVD = 0.1;
    
    if ~exist('enable_full_cpqr','var') || isempty(enable_full_cpqr)
        enable_full_cpqr = 0;
    end
    if ~exist('option','var') || isempty(option)
        option = 'lupp';
    end
    if ~exist('k_max','var') || isempty(k_max)
        k_max = floor(min(m,n)/log2(min(m,n)));
    end
    if k_max >= min(m,n)
        if nargout==2
            [J,VT] = ID_cpqr(A,tol);
        else
            J = ID_cpqr(A,tol);
        end
        return;
    end
    
    if strcmpi(option,'cpqr') % using cpqr with / without sketching
        X = randn(k_max,m)*A;
        [~,R,Jall] = qr(full(X),'vec');
        tmp = diag(R(1:k_max,1:k_max));
        dia = abs(tmp)'.*((1:length(tmp)).^(ETA_L_POW_CPQR));
        k = find(dia<tol,1);
        
        if isempty(k) && enable_full_cpqr % switching to full CPQR (SLOW!)
            if nargout==2
                [J,VT] = ID_cpqr(A,tol);
            else
                J = ID_cpqr(A,tol);
            end
            return;
            
        else
            if isempty(k), k = k_max; end
            J = Jall(1:k);
            if nargout == 2
                VT = A(:,J)\A;
            end
            return;
        end
        
    else % using lupp with binary rank searching
        Xlarge = randn(k_max,m)*A; %(k_max,n)
        s = svd(Xlarge);
        k = find(s'.*((1:length(s)).^(ETA_L_POW_SVD))<tol, 1);
        
        if isempty(k) && enable_full_cpqr % switching to full CPQR (SLOW!)
            if nargout==2
                [J,VT] = ID_cpqr(A,tol);
            else
                J = ID_cpqr(A,tol);
            end
            return;
        else
            if isempty(k), k = k_max; end
            X = Xlarge(1:k,:); % (k,n)
            [~,~,Jall] = lu(full(X'),'vec');
            J = Jall(1:k);
            if nargout == 2
                VT = A(:,J)\A;
            end
            return;
        end
        
    end
end

%%
function [J,VT] = ID_cpqr(A,tol)
    warning('Full CPQR invoked!');
    m = size(A,1);
    n = size(A,2);
    [~,R,Jall] = qr(full(A),'vec');
    tmp = diag(R(1:min(m,n),1:min(m,n)));
    dia = abs(tmp)'.*(1:length(tmp));
    k = find(dia<tol,1);
    if isempty(k), k = length(dia); end

    J = Jall(1:k);
    if nargout == 2
        VT = zeros(k,n);
        R11 = R(1:k,1:k);
        R12 = R(1:k,k+1:end);
        VT(:,Jall) = [eye(k),R11\R12];
    end
end