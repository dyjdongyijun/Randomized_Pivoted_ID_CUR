% clear; close;
% m = 2e4; n = 2e4;
% k = 10; r = 1e3;
% a = 2; b = 1; 
% s = 1e-4;
% rng('default')
% tic;
% target = TargetMatGenerator('snn',m,n,k,r,a,b,s);
% toc
% A = target.A;
% %% 
% save('target.mat','-struct','target')
%%
clear; close;
target = load('target_snn1e45e3-a2b1-k10-r5e3-s2e-4.mat');
A = target.A;
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
time = struct('IDq0',zeros(size(k)),...
              'IDq2',zeros(size(k)),...
              'DEIM',zeros(size(k)),...
              'ACA',zeros(size(k)),...
              'ACAglobal',zeros(size(k)),...
              'CSLUCP',zeros(size(k)),...
              'LS',zeros(size(k)),...
              'OptCUR',zeros(size(k)),...
              'SparseCUR',zeros(size(k)));
%% Optimal CUR
rng('default');
seed = randi(4e9);
eps = 1; 
for t = 1:length(kin)
    rng(seed);
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
rng('default');
%% Sparse Optimal CUR
rng('default');
seed = randi(4e9);
eps = 1; 
for t = 1:length(kin)
    rng(seed);
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
rng('default');
%% 
k = min(min(c),min(r));
randmat = 'gauss'; q = 2; ortho = 1; 

%% CURID q = 0
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'IDq0')
for t = 1:length(k)
    rng(seed);
    tic;
    [i,j] = CUR_ID(A, k(t), randmat, 0);
    time.IDq0(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq0(t) = norm(E); 
    errfro.IDq0(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.IDq0(t))
end
rng('default');

%% CUR-ID q = 2  
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'IDq2')
for t = 1:length(k)
    rng(seed);
    tic;
    [i,j] = CUR_ID(A, k(t), randmat, 2);
    time.IDq2(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.IDq2(t) = norm(E); 
    errfro.IDq2(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.IDq2(t))
end
rng('default');

%% CUR_DEIM
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'CUR-DEIM')
rsvd = 1;
for t = 1:length(k)
    rng(seed);
    tic;
    %     [i,j] = CUR_DEIM(A, k(t));
    [i,j] = CUR_DEIM(A, k(t), rsvd);
    time.DEIM(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.DEIM(t) = norm(E); 
    errfro.DEIM(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.DEIM(t))
end
rng('default');

% CUR_ACA
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'CUR-ACA')
for t = 1:length(k)
    rng(seed);
    tic;
    %     [j,~,i] = CUR_ACA(A, k(t));
    [i,j] = CUR_ACA(A, k(t));
    time.ACA(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.ACA(t) = norm(E); 
    errfro.ACA(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.ACA(t))
end
rng('default');

% CUR_ACA_global
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'CUR-ACA-global')
for t = 1:length(k)
    rng(seed);
    tic;
    %     [j,~,i] = CUR_ACA_global(A, k(t));
    [i,j] = CUR_ACA_global(A, k(t));
    time.ACAglobal(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.ACAglobal(t) = norm(E); 
    errfro.ACAglobal(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.ACAglobal(t))
end
rng('default');

% CUR_LeverageScore
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'CUR-leverage score')
rsvd = 1;
for t = 1:length(k)
    rng(seed);
    tic;
    %     [i,j] = CUR_LeverageScore(A, k(t));
    [i,j] = CUR_LeverageScore(A, k(t), rsvd);
    time.LS(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS(t) = norm(E); 
    errfro.LS(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.LS(t))
end
rng('default');

% CUR_CS_LUCP
rng('default');
seed = randi(4e9);
fprintf('%s \n', 'CUR-CS-LUCP')
for t = 1:length(k)
    rng(seed);
    tic;
    %     [j,~,i] = CUR_CS_LUCP(A, k(t));
    [i,j] = CUR_CS_LUCP(A, k(t));
    time.CSLUCP(t) = toc;
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.CSLUCP(t) = norm(E); 
    errfro.CSLUCP(t) = norm(E,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.CSLUCP(t))
end
rng('default');

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
    sprintf('ID $q=0$: $%.2e$ sec', mean(time.IDq0)),...
    sprintf('ID $q=2$: $%.2e$ sec', mean(time.IDq2)),...
    sprintf('DEIM-RSVD: $%.2e$ sec', mean(time.DEIM)),...
    sprintf('CS-LUCP: $%.2e$ sec', mean(time.CSLUCP)),...
    sprintf('ACA-global: $%.2e$ sec', mean(time.ACAglobal)),...
    sprintf('ACA-os=10: $%.2e$ sec', mean(time.ACA)),...
    sprintf('LS-max: $%.2e$ sec', mean(time.LS)),...
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
%% time plot
semilogy(k, time.IDq0, 'o-', 'LineWidth', 1.5)
hold on
plot(k, time.IDq2, 'o-', 'LineWidth', 1.5)
plot(k, time.DEIM, 's-', 'LineWidth', 1.5)
plot(k, time.CSLUCP, 's-', 'LineWidth', 1.5)
plot(k, time.ACAglobal, 's-', 'LineWidth', 1.5)
plot(k, time.ACA, 'd-', 'LineWidth', 1.5)
plot(k, time.LS, '^-', 'LineWidth', 1.5)
plot(k, time.OptCUR, 'x-', 'LineWidth', 1.5)
plot(k, time.SparseCUR, 'x-', 'LineWidth', 1.5)
hold off
xlabel('$k=min(c,r)$','interpreter','latex')
ylabel('time / sec','interpreter','latex')
legend('ID $q=0$',...
    'ID $q=2$',...
    'DEIM-RSVD',...
    'CS-LUCP',...
    'ACA-global',...
    'ACA-os=10',...
    'LS-max',...
    sprintf('Optimal-CUR-$%.2f$', eps),...
    sprintf('Sparse-Optimal-CUR-$%.2f$', eps),...
    'interpreter','latex')
title(strcat(target.description, '  runtime'),...
    'interpreter','latex')

%% save data
save('rank.mat','k')
save('time.mat','-struct','time')
save('err2.mat','-struct','err2')
save('errfro.mat','-struct','errfro')
