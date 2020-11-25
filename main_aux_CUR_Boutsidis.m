clear; close;
m = 3e4; n = 300;
k = 10; r = n;
a = 1000; b = 1; 
tic;
target = TargetMatGenerator('snn',m,n,k,r,a,b);
toc
A = target.A;
%%
clear; close;
path = 'C:\Users\yijundong\Documents\MATLAB\OdenUT\RandNLA\SparseMatLib\';
filename = strcat(path,'ak2010.mat');
target = load(filename); target = target.Problem;
A = target.A;
%% parameters
kin = 1:10;
% CUR_boutsidis: c and r
c = zeros(size(kin));
r = zeros(size(kin));
%% Output data
err2 = struct('IDq0',zeros(size(k)),...
              'IDq2',zeros(size(k)),...
              'DEIM',zeros(size(k)),...
              'LS',zeros(size(k)),...
              'LSopt',zeros(size(k)));
errfro = struct('IDq0',zeros(size(k)),...
                'IDq2',zeros(size(k)),...
                'DEIM',zeros(size(k)),...
                'LS',zeros(size(k)),...
                'LSopt',zeros(size(k)));
time = struct('IDq0',0,...
              'IDq2',0,...
              'DEIM',0,...
              'LS',0,...
              'LSopt',zeros(size(k)));
%% LS Boutsidis
eps = 1; 
for t = 1:length(kin)
    tic;
    [C,R] = CUR_boutsidis(A,kin(t),eps);
    time.LSopt(t) = toc;
    c(t) = size(C,2);
    r(t) = size(R,1);
    [Qc,~] = qr(C,0); [Qr,~] = qr(R',0);
    CUR = (Qc'*A)*Qr; CUR = (Qc*CUR)*Qr'; 
    E = A - CUR;
    err2.LSopt(t) = norm(full(E)); 
    errfro.LSopt(t) = norm(E,'fro');
    fprintf('progress: %d / %d \n', t, length(kin))
end
fprintf('%s: %.4f\n', 'LSopt', sum(time.LSopt))

%% DEIM + CURID + LS 
k = min(c,r);
randmat = 'gauss'; q = 2; ortho = 1; rsvd = 1;
% CURID q = 0
tic;
[j,~,i] = CUR_ID(A, k(end), randmat, 0);
time.IDq0 = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq0(t) = norm(E); 
    errfro.IDq0(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'IDq0', time.IDq0)
% CURID q = 2
tic;
[j,~,i] = CUR_ID(A, k(end), randmat, 2);
time.IDq2 = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq2(t) = norm(E); 
    errfro.IDq2(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'IDq2', time.IDq2)
% CUR_DEIM
tic;
%     [j,~,i] = CUR_DEIM(A, k(t));
[j,~,i] = CUR_DEIM(A, k(t), rsvd);
time.DEIM = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.DEIM(t) = norm(E); 
    errfro.DEIM(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'DEIM', time.DEIM)
% CUR_LeverageScore
tic;
%     [j,~,i] = CUR_LeverageScore(A, k(t));
[j,~,i] = CUR_LeverageScore(A, k(t), rsvd);
time.LS = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS(t) = norm(E); 
    errfro.LS(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS', time.LS)

%% visualization
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));

subplot(1,2,1)
semilogy(k, sigma(k+1)./sigma(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (err2.IDq0)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.IDq2)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.DEIM)./sigma(1), 's-', 'LineWidth', 1.5)
semilogy(k, (err2.LS)./sigma(1), '^-', 'LineWidth', 1.5)
semilogy(k, (err2.LSopt)./sigma(1), 'x-', 'LineWidth', 1.5)
hold off
xlabel('$k=min(c,r)$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}$',...
    sprintf('ID $q=0$: $%.2e$ sec', time.IDq0),...
    sprintf('ID $q=2$: $%.2e$ sec', time.IDq2),...
    sprintf('DEIM-RSVD: $%.2e$ sec', time.DEIM),...
    sprintf('LS-max: $%.2e$ sec', time.LS),...
    sprintf('LS-Boutsidis-$%.2f$: $%.2e$ sec', eps, mean(time.LSopt)),...
    'interpreter','latex')
title(strcat(target.description, '  spectral norm'),...
    'interpreter','latex')

subplot(1,2,2)
semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (errfro.IDq0)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.IDq2)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.DEIM)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.LS)./sfro(1), '^-', 'LineWidth', 1.5)
semilogy(k, (errfro.LSopt)./sfro(1), 'x-', 'LineWidth', 1.5)
hold off
xlabel('$k=min(c,r)$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}$',...
    'ID $q=0$',...
    'ID $q=2$',...
    'DEIM-RSVD',...
    'LS-max',...
    sprintf('LS-Boutsidis-$%.2f$', eps),...
    'interpreter','latex')
title(strcat(target.description, '  Frobenius norm'),...
    'interpreter','latex')
%% single test
k = 20; eps = 1;
tic;
[C,R] = CUR_sparse_boutsidis(A,k,eps);
sprintf('sparse run time = %.2f sec', toc)
[Qc,~] = qr(C,0);
[Qr,~] = qr(R',0);
sprintf('orthogonalized error = %.2e', normest(A - Qc*((Qc'*A)*Qr)*Qr'))
% sig = target.sigma;
% sig(k+1)
%%
tic;
[C,R] = CUR_boutsidis(B,k,eps);
sprintf('standard run time = %.2f sec', toc)
[Qc,~] = qr(C,0);
[Qr,~] = qr(R',0);
sprintf('orthogonalized error = %.2e', normest(A - Qc*((Qc'*A)*Qr)*Qr'))
%%
