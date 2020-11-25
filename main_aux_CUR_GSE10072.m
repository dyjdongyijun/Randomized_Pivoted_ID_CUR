clear; close;
% load GSE10072 dataset
% gseData = getgeodata('GSE10072', 'ToFile', 'GSE10072.txt');
gseData = geoseriesread('GSE10072.txt');
% struct, fieldnames = {'Header','Data'}
get(gseData.Data)

Araw = double(gseData.Data);
% center data: subtract mean of each row
A = Araw - mean(Araw,2);
% [U,D,V] = svd(A,0);
% A = (U.*(diag(D).^6)')*V';
target = TargetMatGenerator(A,'svd');
A = target.A;
target.description = sprintf('[GSE10072 $%d \\times %d$]',size(A,1),size(A,2));
% save('target.mat','-struct','target')
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
              'LS',zeros(size(k)),...
              'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
errfro = struct('IDq0',zeros(size(k)),...
                'IDq2',zeros(size(k)),...
                'DEIM',zeros(size(k)),...
                'ACA',zeros(size(k)),...
                'LS',zeros(size(k)),...
                'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
time = struct('IDq0',0,...
              'IDq2',0,...
              'DEIM',0,...
              'ACA',0,...
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
%% DEIM + ACA + CURID + LS 
k = min(min(c),min(r));
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
[j,~,i] = CUR_DEIM(A, k(end), rsvd);
time.DEIM = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.DEIM(t) = norm(E); 
    errfro.DEIM(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'DEIM', time.DEIM)

% CUR_ACA
tic;
%     [j,~,i] = CUR_DEIM(A, k(t));
[i,j] = CUR_ACA(A, k(end));
time.ACA = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.ACA(t) = norm(E); 
    errfro.ACA(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'ACA', time.ACA)


% CUR_LeverageScore
tic;
%     [j,~,i] = CUR_LeverageScore(A, k(t));
[j,~,i] = CUR_LeverageScore(A, k(end), rsvd);
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
    sprintf('ACA: $%.2e$ sec', time.ACA),...
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
    'ACA',...
    'LS-max',...
    sprintf('Optimal-CUR-$%.2f$', eps),...
    sprintf('Sparse-Optimal-CUR-$%.2f$', eps),...
    'interpreter','latex')
title(strcat(target.description, '  Frobenius norm'),...
    'interpreter','latex')
%%
