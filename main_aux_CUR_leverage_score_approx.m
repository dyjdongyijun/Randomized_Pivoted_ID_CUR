clear; close;
% load GSE10072 dataset
% gseData = getgeodata('GSE10072', 'ToFile', 'GSE10072.txt');
gseData = geoseriesread('GSE10072.txt');
% struct, fieldnames = {'Header','Data'}
get(gseData.Data)
%%
Araw = double(gseData.Data);
% center data: subtract mean of each row
A = Araw - mean(Araw,2);
% [U,D,V] = svd(A,0);
% A = (U.*(diag(D).^6)')*V';
target = TargetMatGenerator(A,'svd');
A = target.A;
target.description = sprintf('[GSE10072 $%d \\times %d$]',size(A,1),size(A,2));
%% SNN
clear; close;
m = 3e4; n = 300;
k = 10; r = n;
a = 1000; b = 1; 
tic;
target = TargetMatGenerator('snn',m,n,k,r,a,b);
toc
A = target.A;
%% target ranks
k = 1:3:30;
%% CUR-ID, DEIM, LS
randmat = 'gauss'; q = 2; ortho = 1; rsvd = 1;

err2 = struct('LS_presvd',zeros(size(k)),...
              'LS_fast_LS_approximator',zeros(size(k)),...
              'LS_RSVD',zeros(size(k)),...
              'LS_SVD',zeros(size(k)));
errfro = struct('LS_presvd',zeros(size(k)),...
              'LS_fast_LS_approximator',zeros(size(k)),...
              'LS_RSVD',zeros(size(k)),...
              'LS_SVD',zeros(size(k)));
time = struct('LS_presvd',0,...
              'LS_fast_LS_approximator',0,...
              'LS_RSVD',0,...
              'LS_SVD',0);


% LS_presvd
tic;
[j,~,i] = CUR_LeverageScore(A, k(end), 0, target.U, target.V);
time.LS_presvd = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS_presvd(t) = norm(E); 
    errfro.LS_presvd(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS_presvd', time.LS_presvd)

% LS_fast_LS_approximator
tic;
[j,~,i] = CUR_LeverageScore(A, k(end), 2);
time.LS_fast_LS_approximator = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS_fast_LS_approximator(t) = norm(E); 
    errfro.LS_fast_LS_approximator(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS_fast_LS_approximator', time.LS_fast_LS_approximator)

% LS_RSVD
tic;
[j,~,i] = CUR_LeverageScore(A, k(end), 1);
time.LS_RSVD = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS_RSVD(t) = norm(E); 
    errfro.LS_RSVD(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS_RSVD', time.LS_RSVD)

% LS_SVD
tic;
[j,~,i] = CUR_LeverageScore(A, k(end), 0);
time.LS_SVD = toc;
for t = 1:length(k)
    E = CUR_Error(A,i(1:k(t)),j(1:k(t)));
    err2.LS_SVD(t) = norm(E); 
    errfro.LS_SVD(t) = norm(E,'fro');
end
fprintf('%s: %.4f\n', 'LS_SVD', time.LS_SVD)

%%
% save('smprk_2.mat', 'k');
% save('target_2_snn3e4x300a1000b1.mat', '-struct', 'target');
% save('err2_2_snn3e4x300a1000b1.mat', '-struct', 'err2');
% save('errfro_2_snn3e4x300a1000b1.mat', '-struct', 'errfro');
% save('time_2_.mat', '-struct', 'time');
% %% 
% k = load('smprk_2.mat'); k = k.k;
% target = load(''); 
% err2 = load('err2_2_.mat');
% errfro = load('errfro_2_.mat');
%% visualization
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));

subplot(1,2,1)
semilogy(k, sigma(k+1)./sigma(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (err2.LS_presvd)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.LS_fast_LS_approximator)./sigma(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (err2.LS_RSVD)./sigma(1), 's-', 'LineWidth', 1.5)
semilogy(k, (err2.LS_SVD)./sigma(1), 's-', 'LineWidth', 1.5)
hold off
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}$',...
    sprintf('LS presvd: $%.2e$ sec', time.LS_presvd),...
    sprintf('LS fast LS approximator: $%.2e$ sec', time.LS_fast_LS_approximator),...
    sprintf('LS RSVD: $%.2e$ sec', time.LS_RSVD),...
    sprintf('LS SVD: $%.2e$ sec', time.LS_SVD),...
    'interpreter','latex')
title(strcat(target.description, '  spectral norm'),...
    'interpreter','latex')

subplot(1,2,2)
semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (errfro.LS_presvd)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.LS_fast_LS_approximator)./sfro(1), 'o-', 'LineWidth', 1.5)
semilogy(k, (errfro.LS_RSVD)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.LS_SVD)./sfro(1), 's-', 'LineWidth', 1.5)
hold off
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}$',...
    'LS presvd',...
    'LS fast LS approximator',...
    'LS RSVD',...
    'LS SVD',...
    'interpreter','latex')
title(strcat(target.description, '  Frobenius norm'),...
    'interpreter','latex')
