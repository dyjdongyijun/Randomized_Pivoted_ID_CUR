function [I,J,U] = CUR_LocalSketchACA(A,k)
%   A = [m*n matrix]
%   k = target rank
%   I= [1*k subset of [m]] row skeletons
%   J= [1*k subset of [n]] column skeletons
%   U = [k*k] matrix
    os = 5; 
    maxiter = 10;
    p = k + os; 
    q = p + os;
%   row skeleton
    [Iout, kout] = rowID_LocalSketchACA(A,q,p);
    i = 0;
    while kout < p && i < maxiter
        [Iout, kout] = rowID_LocalSketchACA(A,q,p);
        i = i+1;
    end
    I = Iout(1:p);
%   column skeleton
    [Jout, kout] = rowID_LocalSketchACA(A(I,:)',p,k);
    i = 0;
    while kout < k && i < maxiter
        [Jout, kout] = rowID_LocalSketchACA(A(I,:)',p,k);
        i = i+1;
    end
    J = Jout(1:k);
    
    I = I(1:k); J = J(1:k);
    if nargout==3
        C = A(:,J(1:k)); R = A(I(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end
%%
function [I, kout] = rowID_LocalSketchACA(A,q,k)
%   A = [m*n]
%   q = embedding dim (Y = A*Pi = [m*q])
%   k = target rank
%   I = [1*kout subset of [m]] row skeletons
    eps = 1e-10;

    [m,n] = size(A);
    U = zeros(m,k); V = zeros(n,k);
    I = 1:m; J = 1:q;
    
    Pi = randn(n,q);
    j = randi(q); 
    Iact = I; Jact = J;
    
    for t = 1:k
        if ~(j==1)
            J(j+t-1) = J(t);
            J(t) = Jact(j);
        end
        
        u = A(Iact,:)*Pi(:,Jact(j)) - U(Iact,1:t-1)*V(Jact(j),1:t-1)';
        [umax,i] = max(abs(u)); 
        if umax < eps 
            t = t-1;
            break; 
        end
        if ~(i==1)
            I(i+t-1) = I(t);
            I(t) = Iact(i);
        end
        u = u./u(i);
        v = (A(Iact(i),:)*Pi(:,Jact))' - V(Jact,1:t-1)*U(Iact(i),1:t-1)';
        
%         dfro = sum(u.^2)*sum(v.^2);
%         fro = fro + dfro + sum((V(Jact,1:t-1)'*v).*(U(Iact,1:t-1)'*u));  
        
        U(Iact,t) = u; 
        V(Jact,t) = v;
        
        Iact = I(t+1:end);
        Jact = J(t+1:end);
        [~,j] = max(abs(V(Jact,t))); 
    end
    
    kout = t;
    I = I(1:k); 
%     J = J(1:k); U = U(:,1:kout); V = V(:,1:kout);
end
