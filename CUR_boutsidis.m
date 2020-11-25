function [C,R,U] = CUR_boutsidis(A,k,eps)
% % %
% implementation of Algorithm 2 in Boutsidis Woodruff 2014 
% arxiv: 1405.7910
% INPUT:
% A, m by n matrix
% k, rank parameter, smaller than rank(A)
% eps, accuracy parameter between 0 and 1
% OUTPUT: CUR decomposition of A
% C = [m*c], R = [r*n]; c \neq r; c,r = O(k + k/eps) 
% U = [c*r]
% % %

%% step1
p = 1;
[~,~,Z1] = rsvd(A,k,p);
h1 = ceil(16*k*log(20*k));

[indOmega1,D1] = randSampling(Z1,h1);
% M1 = Z1'*Omega1*D1;
M1 = Z1';
M1 = M1(:,indOmega1) .* D1';
[~,~,v] = svd(M1,'econ');

% temp = (A - A*(Z1*Z1'))*Omega1*D1;
temp = A - A*Z1*Z1';
temp = temp(:,indOmega1) .* D1';

[s,ind1] = bssSampling(v,temp',4*k);
% C1 = A*Omega1*D1;
C1 = A(:,indOmega1) .* D1';
C1 = C1(:,ind1) .* s';

% c2 = ceil(1620*k/eps); % constant too big here
c2 = ceil(k/eps);
[C2,ind2] = adaptiveCols(A,C1,c2);
C = [C1 C2];

%% step 2
[Y,Delta] = bestSubspaceSVD(A,C,k);
B = Y*Delta;
[Z2,D] = qr(B,0);
E2 = A' - A'*Z2*Z2';

h2 = ceil(8*k*log(20*k));
[indOmega2,D2] = randSampling(Z2,h2);
% M2 = Z2'*Omega*D2;
M2 = Z2';
M2 = M2(:,indOmega2) .* D2';
[~,~,v2] = svd(M2,'econ');

temp = E2(:,indOmega2) .* D2';
[s2,inds2] = bssSampling(v2,temp',4*k);
R1 = A';
R1 = R1(:,indOmega2) .* D2';
R1 = R1(:,inds2) .* s2';
R1 = R1';

r2 = ceil(k/eps);
R2 = adaptiveRows(A,R1,r2);
R = [R1' R2']';

%% step 3
    if nargout==3
%         C = A(:,J(1:k)); R = A(I(1:k),:);
        U = CUR_U_eval(A,C,R);
    end
end

function [Y,Delta] = bestSubspaceSVD(A,V,k)
    [Y,~] = qr(V,0);
    Theta = full(Y'*A);
    [u,~,~] = svd(Theta,'econ');
    Delta = u(:,1:k);
end


function [s,ind] = bssSampling(V,Acal,r)
% % %
% Algorithm 1 in NEAR-OPTIMAL COLUMN-BASED MATRIX RECONSTRUCTION
% with modification in section 8
% INPUT: V = [v1,...,vn], VV' = Id_k
% % %
V = V';
Acal = Acal';
[k,n] = size(V);
s = zeros(n,1);
A = zeros(k,k);

delta_L = 1;
delta_U = (norm(Acal,'fro'))^2/(1-sqrt(k/r));
for tau = 0:r-1
   L = tau - sqrt(r*k);
   % U = delta_U*tau;
   j = 1;
   while j <= n
       U_F = fcnU_F(Acal(:,j),delta_U);
       ll = fcnL(V(:,j),delta_L,A,L);
       if (U_F <= ll)||(j==n)
           break;
       else
           j = j + 1;
       end
   end
   t = 2/(U_F+ll);
   s(j) = s(j) + t;
   temp = V(:,j)*(V(:,j))';
   A = A + t*temp;
end
s = 1/r*(1-sqrt(k/r))*s;

% construct S
s = sqrt(s);
ind = 1:n;
ind = ind(s~=0);
s = s(s~=0);
% S = sparse(ind,1:length(ind),s,n,length(ind));

end


function y = fcnL(v,delta_L,A,L)
temp = (A - (L+delta_L)*eye(size(A))) \ v;
num = (A - (L+delta_L)*eye(size(A))) \ temp;
num = v'*num;

eigA = eig(A);
denom = sum(1./(eigA-L-delta_L) - 1./(eigA-L) );
y = num/denom - v'*temp;
end

function y = fcnU_F(a,delta_U)
y = 1/delta_U*(a'*a);
end

function [rind,D] = randSampling(X,r)
% % %
% this is in definition 6 of Clarkson, Woodruff 2013 paper
% % %
p = sum(X.*X,2);
p = sqrt(p);
p = p/sum(p);
[n,~] = size(X);
rind = zeros(r,1);
D = zeros(r,1);
for i = 1:r
    rind(i) = rsample(1:n,p);
    D(i) = 1/sqrt(p(rind(i))*r);
end
end

function [R2,ind] = adaptiveRows(A,R1,r2)
% % %
% lemma 3.10 of arxiv 1405.7910
% % % 
[~,n] = size(A);
B = A - A*(R1 \ R1);
p = sum(B.^2);
p = p/sum(p);
ind = zeros(r2,1);
for i = 1:r2
    ind(i) = rsample(1:n,p);
end
R2 = A(ind,:);
end

function [C2,ind] = adaptiveCols(A,V,c2)
% % %
% lemma 3.9 of arxiv 1405.7910
% % % 
[~,n] = size(A);
B = A - V*(V\A);
p = sum(B.^2);
p = p/sum(p);
ind = zeros(c2,1);
for i = 1:c2
    ind(i) = rsample(1:n,p);
end
C2 = A(:,ind);
end

function [U,S,V] = rsvd(A,k,p)
% % %
% randomized svd decomposition of A
% % %
% Y = srht(A,k+p);
Y = SRFT(A,k+p);
[Q,~] = qr(Y,0);

B = Q'*A;
[U,S,V] = svd(B,'econ');
U = U(:,1:k);
S = S(1:k,1:k);
V = V(:,1:k);
U = Q*U;
end

function y = rsample(x,p)
cum_p = cumsum(p);
r = rand;
y = x(find(r<cum_p,1,'first'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code below are not working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,D,W] = sparseSVD(A,k,eps)
% % %
% Theorem 47 in Clarkson, Woodruff 2013 paper
% "Low Rank Approximation and Regression in Input Sparsity Time"
% A, n by d matrix
% % %
%% step 1
t = ceil(k^2*(log(k/eps))^6 + k/eps);
Y = sparseEmbed(A,t);
tt = ceil(k/eps*log(k/eps));
Y = srht(Y,tt);
[U,~] = qr(Y,0);
%% step 2
v = ceil(k^2*eps^(-4)*(log(k/eps))^6);
vv = ceil(k*eps^(-3)*(log(k/eps))^2);
SU = sparseEmbed(U',v);
SU = srht(SU,vv);
SU = SU';
SA = sparseEmbed(A',v);
SA = srht(SA,vv);
SA = SA';
%% step 3
[Ut,Sigt,Vt] = svd(SU,'econ');
%% step 4
temp = Ut' * SA;
[u,s,v] = svd(temp,'econ');
temp = u(:,1:k) * s(1:k,1:k) * (v(:,1:k))';
temp = Vt*(Sigt)*temp;
[Uhat,D,W] = svd(temp,'econ');
L = U*Uhat;

end

function C = sparseEmbed(A,t)
% % %
% % Theorem 47 in Clarkson, Woodruff 2013 paper
% "Low Rank Approximation and Regression in Input Sparsity Time"
% A, n by d matrix
% % %
n = size(A,2);
sgn = randi(2,[1,n])*2 - 3; % sampled +1 and -1
A = bsxfun(@times,A,sgn); % flip signs of column of A
idx = sort(randsample(1:n,t,true));
C = A(:,idx); % subsampling
end

