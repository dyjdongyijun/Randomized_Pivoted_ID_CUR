function E = ID_Error(A,C,V)
    if exist('V','var')
        E = A - C*V';
    else
        if issparse(A) 
            [B, T] = qr(C, A, 0); % B:(k,n), T:(k,k), T\B = least_sq(C\A):(k,n)
            VT = T\B;
        else
            VT = C\A; % X = [k*n]
        end
        E = A - C*VT;
    end
end