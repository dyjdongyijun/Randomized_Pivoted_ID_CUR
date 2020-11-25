function [C,R,U] = CUR_sparse_boutsidis(A,k,eps)
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
% [~,~,Z1] = sparseSVD(A,k,eps); %Z1 = [n*k]
[~,~,Z1] = RSVD(A,k,'sparse'); %Z1 = [n*k]
h1 = ceil(16*k*log(20*k)); % h1>k

[indOmega1,D1] = randSampling(Z1,h1);
% M1 = Z1'*Omega1*D1;
M1 = Z1(indOmega1,:)' .* D1'; % M1 = [k*h1]
[~,~,v] = svd(M1,'econ'); % v = [h1*k]

% temp = (A - A*(Z1*Z1'))*Omega1*D1;
temp = A - A*Z1*Z1'; % [m*n]
temp = temp(:,indOmega1) .* D1'; % [m*h1]

% [s,ind1] = bssSampling(v,temp',4*k); % v = [h1*k], temp' = [h1*m], sample 4k out of h1
[s,ind1] = bssSamplingSparse(v,temp',4*k,1);
% C1 = A*Omega1*D1;
C1 = A(:,indOmega1) .* D1';
C1 = C1(:,ind1) .* s';

% c2 = ceil(1620*k/eps); % constant too big here
c2 = ceil(k/eps);
[C2,ind2] = adaptiveColsSparse(A,C1,c2);
C = [C1 C2];

%% step 2
[Y,Delta] = approxSubspaceSVD(A,C,k);
B = Y*Delta;
[Z2,D] = qr(B,0);
E2 = A' - (A'*Z2)*Z2';

h2 = ceil(8*k*log(20*k));
[indOmega2,D2] = randSampling(Z2,h2);
% M2 = Z2'*Omega*D2;
M2 = Z2';
M2 = M2(:,indOmega2) .* D2';
[~,~,v2] = svd(M2,'econ');

temp = E2(:,indOmega2) .* D2';
% [s2,inds2] = bssSampling(v2,temp',4*k);
[s2,inds2] = bssSamplingSparse(v2,temp',4*k,1);
R1 = A';
R1 = R1(:,indOmega2) .* D2';
R1 = R1(:,inds2) .* s2';
R1 = R1';

r2 = ceil(k/eps);
R2 = adaptiveRowsSparse(A,R1,r2);
R = [R1' R2']';
%%
if nargout==3
    U = CUR_U_eval(A,C,R);
end
end

%%
function [s,ind] = bssSamplingSparse(V,A,r,eps)
% % %
% Algorithm 1 in NEAR-OPTIMAL COLUMN-BASED MATRIX RECONSTRUCTION
% with modification in section 8
% INPUT: V = [v1,...,vn], VV' = Id_k
% % % 
    [n,l] = size(A);
    xi = min((n/eps)^2, l);
    B = sparseEmbed(A', xi)';
    [s,ind] = bssSampling(V,B,r);
end

%%
function [s,inz] = bssSampling(V,Acal,r)
% % %
% Algorithm 1 in NEAR-OPTIMAL COLUMN-BASED MATRIX RECONSTRUCTION
% with modification in section 8
% INPUT: V = [v1;...;vn], V'*V = Id_k, V = [n*k]
%        Acal = [a1;...;an] = [n*l]
%        k < r < n
% OUTPUT: s = [r*1] weights of nonzeros
%         inz = [r*1] index of nonzeros
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
inz= (1:n)';
inz = inz(s~=0);
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
%%
function [rind,D] = randSampling(X,r)
% % %
% this is in definition 6 of Clarkson, Woodruff 2013 paper
% rind = [r*1] s.t. Omega = X(rind,:)
% D = [r*1] diagonal of the rescaling matrix
% % %
p = sum(X.^2,2);
% p = sqrt(p);
p = p./sum(p);
n = size(X,1);
rind = sampleProb((1:n)', p, r);
D = 1./sqrt(p(rind).*r);
end

%%
function [R2,ind] = adaptiveRowsSparse(A,R1,r2)
% % %
% lemma 3.10 of arxiv 1405.7910
% A = [m*n], R1 = [r1*n], R2 = [r2*n] = A(ind,:)
% % % 
[m,~] = size(A);
Baux = A - A*(R1 \ R1); % Baux = [m*n]
B = JLT(Baux', 1)'; % B = [m*s], s<<n
p = sum(B.^2,2); % p = [m*1]
p = p/sum(p);
ind = sampleProb((1:m)', p, min(r2,m));
R2 = A(ind,:);
end

%%
function [C2,ind] = adaptiveColsSparse(A,V,c2)
% % %
% lemma 3.9 of arxiv 1405.7910
% A = [m*n], V = [m*c1], C2 = [m*c2] = A(:,ind)
% % % 
[~,n] = size(A);
Baux = A - V*(V\A); % Baux = [m*n]
B = JLT(Baux, 1); % B = [s*n], s<<m
p = sum(B.^2);
p = p/sum(p);
ind = sampleProb((1:n)', p, min(c2,n));
C2 = A(:,ind);
end

%%
function Bt = JLT(B, beta)
%   B = [m*n], S = [s*m]
%   Bt = SB = s*n
    [m,n] = size(B);
    s = floor((4 + 2*beta)/(1/2^2-1/2^3)*log(n));
    S = (-1).^randi(2,[s,m])./sqrt(s);
    Bt = S*B;
end

%%
function [Y,Delta] = bestSubspaceSVD(A,V,k)
    [Y,~] = qr(V,0);
    Theta = full(Y'*A);
    [u,~,~] = svd(Theta,'econ');
    Delta = u(:,1:k);
end

%%
function [Y,Delta] = approxSubspaceSVD(A,V,k,eps)
%   A = [m*n], V = [m*c]
%   Y = [m*c], Delta = [c*k]
    if nargin<4, eps = 1; end
    [Y,~] = qr(V,0);
    xi = size(V,2)^2/eps^2;
    Xi = sparseEmbed(A'*Y, min(xi, size(A,2)))'; % Xi = [c*xi]
    [u,~,~] = svd(full(Xi),'econ');
    Delta = u(:,1:k);
end

%%
% % %
% sparseSVD does not perform as well as RSVD (with gaussian random
% matrices) in practice
% % % 
function [L,D,W] = sparseSVD(A,k,eps)
% % %
% Theorem 47 in Clarkson, Woodruff 2013 paper
% "Low Rank Approximation and Regression in Input Sparsity Time"
% A = [n*d]
% % %
%% step 1
t = min(ceil(k^2*(log(k/eps))^6 + k/eps), 10*k);
tt = min(ceil(k/eps*log(k/eps)), t);
Y = sparseEmbed(A',t)';
Y = SRFT(Y,tt);
[U,~] = qr(Y,0);
%% step 2
v = min(ceil(k^2*eps^(-4)*(log(k/eps))^6), 10*k);
vv = min(ceil(k*eps^(-3)*(log(k/eps))^2), v);
SU = sparseEmbed(U,v);
SU = SRFT(SU',vv)';
SA = sparseEmbed(A,v);
SA = SRFT(SA',vv)';
%% step 3
[Ut,Sigt,Vt] = svd(SU,'econ'); sigt = diag(Sigt);
%% step 4
[u,s,v] = svd(Ut'*SA,'econ'); s = diag(s); 
aux = u(:,1:k)*(s(1:k).*v(:,1:k)');
aux = Vt*((1./sigt).*aux);
[Uh,D,W] = svd(aux,'econ');
L = U*Uh;
end


%%
function R = sparseEmbed(A,t)
% % %
% % Theorem 47 in Clarkson, Woodruff 2013 paper
% "Low Rank Approximation and Regression in Input Sparsity Time"
% Subspace embedding of Column space n->t
% A = [n*d]
% R = (Phi*D)*A = [t*d]
% Phi = [t*n]
% % %
n = size(A,1);
d = (-1).^randi(2,[n,1]); % sampled +1 and -1
DA = d.*A; % flip signs of column of A -> O(n*d)
h = randi(t,[1,n]); % hash function
Phi = sparse(h,1:n,ones(1,n));
R = Phi*DA; % subsampling
if issparse(A), R = sparse(R); end
end
