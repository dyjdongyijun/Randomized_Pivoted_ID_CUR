% random sparse matrix
clear; close;
m = 3e4; n = 1e3;
k = 10; r = n;
a = 2; b = 1; 
tic;
target = TargetMatGenerator('snn',m,n,k,r,a,b);
toc
A = target.A;
%%
% % %
% All test matrices are from the UF sparse matrix collections
% bips07_1998.mat
% % % 
clear; close;
path = 'C:\Users\yijundong\Documents\MATLAB\OdenUT\RandNLA\SparseMatLib\';
filename = strcat(path,'bips07_1998.mat');
target = load(filename); target = target.Problem;
A = target.A;
target.description = sprintf('bips07_1998: %d * %d', size(A,1), size(A,2));
target.sigma = svds(A,51);

%% single test
k = 20; eps = 1;
tic;
[C,R] = CUR_sparse_boutsidis(A,k,eps);
sprintf('sparse run time = %.2f sec', toc)
% bips07_1998.mat: sparse run time = 33.08 sec
[Qc,~] = qr(C,0);
[Qr,~] = qr(R',0);
sprintf('orthogonalized error = %.2e', normest(A - Qc*((Qc'*A)*Qr)*Qr'))
% bips07_1998.mat: orthogonalized error = 5.80e+04

tic;
[C,R] = CUR_boutsidis(A,k,eps);
sprintf('standard run time = %.2f sec', toc) 
% bips07_1998.mat: standard run time = 81.41 sec
[Qc,~] = qr(C,0);
[Qr,~] = qr(R',0);
sprintf('orthogonalized error = %.2e', normest(A - Qc*((Qc'*A)*Qr)*Qr'))
% bips07_1998.mat: orthogonalized error = 5.20e+04

inf = min(size(C,2), size(R,1));
sup = max(size(C,2), size(R,1));
sig = svds(A, sup+1);
sprintf('best rank-k error = %.2e', sig(inf+1))
% bips07_1998.mat: best rank-k error = 2.50e+04

%% parameters
kin = 1:10;
% CUR_boutsidis: c and r
c = zeros(2,length(kin)); % [OptCUR; SparseCUR]
r = zeros(2,length(kin));
%% Output data
k = kin;
err2 = struct('IDq0',zeros(size(k)),...
              'IDq2',zeros(size(k)),...
              'DEIM',zeros(size(k)),...
              'ACA',zeros(size(k)),...
              'ACAglobal',zeros(size(k)),...
              'CSLUCP',zeros(size(k)),...
              'LS',zeros(size(k)),...
              'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
errfro = struct('IDq0',zeros(size(k)),...
                'IDq2',zeros(size(k)),...
                'DEIM',zeros(size(k)),...
                'ACA',zeros(size(k)),...
                'ACAglobal',zeros(size(k)),...
                'CSLUCP',zeros(size(k)),...
                'LS',zeros(size(k)),...
                'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
time = struct('IDq0',0,...
              'IDq2',0,...
              'DEIM',0,...
              'ACA',0,...
              'ACAglobal',0,...
              'CSLUCP',0,...
              'LS',0,...
              'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
%% Optimal CUR
eps = 1; 
for t = 1:length(kin)
    tic;
    [C,R] = CUR_boutsidis(A,kin(t),eps);
    time.OptCUR(t) = toc;
    c(1,t) = size(C,2);
    r(1,t) = size(R,1);
    [Qc,~] = qr(C,0); [Qr,~] = qr(R',0);
    CUR = (Qc'*A)*Qr; CUR = (Qc*CUR)*Qr'; 
    E = A - CUR;
    err2.OptCUR(t) = norm(full(E)); 
    errfro.OptCUR(t) = norm(E,'fro');
    fprintf('progress: %d / %d -> (%.2e sec) \n', t, length(kin), time.OptCUR(t))
end
%% Sparse Optimal CUR
eps = 1; 
for t = 1:length(kin)
    tic;
    [C,R] = CUR_sparse_boutsidis(A,kin(t),eps);
    time.SparseCUR(t) = toc;
    c(2,t) = size(C,2);
    r(2,t) = size(R,1);
    [Qc,~] = qr(C,0); [Qr,~] = qr(R',0);
    CUR = (Qc'*A)*Qr; CUR = (Qc*CUR)*Qr'; 
    E = A - CUR;
    err2.SparseCUR(t) = norm(full(E)); 
    errfro.SparseCUR(t) = norm(E,'fro');
    fprintf('progress: %d / %d -> (%.2e sec) \n', t, length(kin), time.SparseCUR(t))
end
%% CUR-ID q = 0
k = min(min(c),min(r));
randmat = 'gauss'; q = 2; ortho = 1; 
% CURID q = 0
tic;
[i,j] = CUR_ID(A, k(end), randmat, 0);
time.IDq0 = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq0(t) = norm(E); 
    errfro.IDq0(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'IDq0', time.IDq0)

% CUR-ID q = 2   
tic;
[i,j] = CUR_ID(A, k(end), randmat, 2);
time.IDq2 = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq2(t) = norm(E); 
    errfro.IDq2(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'IDq2', time.IDq2)

% CUR_DEIM
rsvd = 1;
tic;
%     [i,j] = CUR_DEIM(A, k(t));
[i,j] = CUR_DEIM(A, k(end), rsvd);
time.DEIM = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.DEIM(t) = norm(E); 
    errfro.DEIM(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'DEIM', time.DEIM)

% CUR_ACA
tic;
%     [j,~,i] = CUR_ACA(A, k(t));
[i,j] = CUR_ACA(A, k(end));
time.ACA = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.ACA(t) = norm(E); 
    errfro.ACA(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'ACA', time.ACA)

% CUR_ACA_global

% err2.ACAglobal = zeros(size(k));
% errfro.ACAglobal = zeros(size(k));
% time.ACAglobal = 0;

tic;
%     [j,~,i] = CUR_ACA_global(A, k(t));
[i,j] = CUR_ACA_global(A, k(end));
time.ACAglobal = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.ACAglobal(t) = norm(E); 
    errfro.ACAglobal(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'ACAglobal', time.ACAglobal)

% CUR_LeverageScore
rsvd = 1;
tic;
%     [i,j] = CUR_LeverageScore(A, k(t));
[i,j] = CUR_LeverageScore(A, k(end), rsvd);
time.LS = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS(t) = norm(E); 
    errfro.LS(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS', time.LS)

%% CUR_CS_LUCP

% err2.CSLUCP = zeros(size(k));
% errfro.CSLUCP = zeros(size(k));
% time.CSLUCP = 0;

tic;
%     [j,~,i] = CUR_CS_LUCP(A, k(t));
[i,j] = CUR_CS_LUCP(A, k(end));
time.CSLUCP = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.CSLUCP(t) = norm(E); 
    errfro.CSLUCP(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'CSLUCP', time.CSLUCP)

%% visualization
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));

subplot(1,2,1)
semilogy(k, sigma(k+1)./sigma(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (err2.IDq0)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.IDq2)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.DEIM)./sigma(1), 's-', 'LineWidth', 1.5)
semilogy(k, (err2.CSLUCP)./sigma(1), 's-', 'LineWidth', 1.5)
semilogy(k, (err2.ACAglobal)./sigma(1), 's-', 'LineWidth', 1.5)
semilogy(k, (err2.ACA)./sigma(1), 'd-', 'LineWidth', 1.5)
semilogy(k, (err2.LS)./sigma(1), '^-', 'LineWidth', 1.5)
semilogy(k, (err2.OptCUR)./sigma(1), 'x-', 'LineWidth', 1.5)
semilogy(k, (err2.SparseCUR)./sigma(1), 'x-', 'LineWidth', 1.5)
hold off
xlabel('$k=min(c,r)$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_1$',...
    sprintf('ID $q=0$: $%.2e$ sec', time.IDq0),...
    sprintf('ID $q=2$: $%.2e$ sec', time.IDq2),...
    sprintf('DEIM-RSVD: $%.2e$ sec', time.DEIM),...
    sprintf('CS-LUCP: $%.2e$ sec', time.CSLUCP),...
    sprintf('ACA-global: $%.2e$ sec', time.ACAglobal),...
    sprintf('ACA-os=10: $%.2e$ sec', time.ACA),...
    sprintf('LS-max: $%.2e$ sec', time.LS),...
    sprintf('Optimal-CUR-$%.2f$: $%.2e$ sec', eps, mean(time.OptCUR)),...
    sprintf('Sparse-Optimal-CUR-$%.2f$: $%.2e$ sec', eps, mean(time.SparseCUR)),...
    'interpreter','latex')
title(strcat(target.description, '  spectral norm'),...
    'interpreter','latex')

subplot(1,2,2)
semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (errfro.IDq0)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.IDq2)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.DEIM)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSLUCP)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.ACAglobal)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.ACA)./sfro(1), 'd-', 'LineWidth', 1.5)
semilogy(k, (errfro.LS)./sfro(1), '^-', 'LineWidth', 1.5)
semilogy(k, (errfro.OptCUR)./sfro(1), 'x-', 'LineWidth', 1.5)
semilogy(k, (errfro.SparseCUR)./sfro(1), 'x-', 'LineWidth', 1.5)
hold off
xlabel('$k=min(c,r)$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sigma_1$',...
    'ID $q=0$',...
    'ID $q=2$',...
    'DEIM-RSVD',...
    'CS-LUCP',...
    'ACA-global',...
    'ACA-os=10',...
    'LS-max',...
    sprintf('Optimal-CUR-$%.2f$', eps),...
    sprintf('Sparse-Optimal-CUR-$%.2f$', eps),...
    'interpreter','latex')
title(strcat(target.description, '  Frobenius norm'),...
    'interpreter','latex')
%%
