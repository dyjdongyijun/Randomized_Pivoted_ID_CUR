% larger test
% clear; close;
% tag = 'p2p-Gnutella09';
% path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
% path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
% ranks = load(strcat(path, sprintf('rank_%s.mat',tag))); ranks = ranks.k;
% target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\lapack')
%% {DetCPQR, LUCP} for {'full', 'svd', 'sketch'} for column ID
clear; close;
% seed = {'gauss', 500, 500, 500, @(x) log(x)}; % {'gauss', m, n, k, f}
% seed = {'gauss', 500, 500, 500, @(x) x.^2}; % {'gauss', m, n, k, f}
% seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
seed = {'laplacian', floor(500), floor(1000)}; % {'laplacian',n,m}
target = TargetMatGenerator(seed{:});
% tag = 'small-gaussian-m500-n500-r500-logdecay';
% tag = 'small-gaussian-m500-n500-r500-polydecay';
tag = 'small-laplacian-n500-m1000';
plot(target.sigma, 'k.-')
ranks = 20:20:400;
%%
A = target.A;
fprintf('%s \n', tag)
cmp = {'full', 'svd', 'sketch'};

err2_cpqr = struct();
errfro_cpqr = struct();
err2_lucp = struct();
errfro_lucp = struct();
for i = 1:length(cmp)
    err2_cpqr.(cmp{i}) = zeros(size(ranks));
    errfro_cpqr.(cmp{i}) = zeros(size(ranks));
    err2_lucp.(cmp{i}) = zeros(size(ranks));
    errfro_lucp.(cmp{i}) = zeros(size(ranks));
end

%% CPQR
% full: DetCPQR
base = 'full'; 
k = ranks(end);
tic;
[j,~] = DetColumnID(A,k);
[i, ~] = DetColumnID(A(:,j(1:k))',k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_cpqr.(base)(idx) = normest(E);
    else
        err2_cpqr.(base)(idx) = norm(E);
    end
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% svd
base = 'svd';
U = target.U; % (m,n), r=min(m,n)
W = target.V; % (n,r), r=min(m,n)
k = ranks(end);
tic;
[j,~] = DetColumnID(W(:,1:k)',k);
[i,~] = DetColumnID(U(:,1:k)',k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_cpqr.(base)(idx) = normest(E);
    else
        err2_cpqr.(base)(idx) = norm(E);
    end
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% sketch: Stream-CPQR
base = 'sketch';
os = 10;
m = size(A,1); n = size(A,2);
k = ranks(end);
tic;
Sm = embed(m, min(k+os,m));
X = Sm(A); % (k+os, n)
Sn = embed(n, min(k+os,n));
YT = Sn(A'); % (k+os, m)
[j,~] = DetColumnID(X,k);
[i,~] = DetColumnID(YT,k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_cpqr.(base)(idx) = normest(E);
    else
        err2_cpqr.(base)(idx) = norm(E);
    end
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

err2file = sprintf('err2-base-cpqr_%s', tag);
errfrofile = sprintf('errfro-base-cpqr_%s', tag);
save(err2file,'-struct','err2_cpqr')
save(errfrofile,'-struct','errfro_cpqr')
%% LUCP
% full: LUCP
base = 'full';
k = ranks(end);
tic;
[i,j] = CUR_LUCP(full(A),k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_lucp.(base)(idx) = normest(E);
    else
        err2_lucp.(base)(idx) = norm(E);
    end
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% svd
base = 'svd';
U = target.U; % (m,r), r=min(m,n)
W = target.V; % (n,r), r=min(m,n)
k = ranks(end);
tic;
[j,~] = CUR_LUCP(W(:,1:k),k);
[i,~] = CUR_LUCP(U(:,1:k),k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_lucp.(base)(idx) = normest(E);
    else
        err2_lucp.(base)(idx) = norm(E);
    end
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% sketch
base = 'sketch';
os = 10;
m = size(A,1);
k = ranks(end);
tic;
Sm = embed(m, min(k+os,m));
X = Sm(A); % (k+os, n)
Sn = embed(n, min(k+os,n));
Y = Sn(A')'; % (m, k+os)
[j,~] = CUR_LUCP(full(X'),k);
[i,~] = CUR_LUCP(full(Y),k);
t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    C = A(:,j(1:k));
    R = A(i(1:k),:);
    E = CUR_Error(A,C,R);
    if issparse(E)
        err2_lucp.(base)(idx) = normest(E);
    else
        err2_lucp.(base)(idx) = norm(E);
    end
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% err2file = sprintf('err2-base-lucp_%s', tag);
% errfrofile = sprintf('errfro-base-lucp_%s', tag);
% save(err2file,'-struct','err2_lucp')
% save(errfrofile,'-struct','errfro_lucp')
%% load
path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
% err2_cpqr = load(strcat(path, sprintf('err2-base-cpqr_%s.mat',tag)));
% errfro_cpqr = load(strcat(path, sprintf('errfro-base-cpqr_%s.mat',tag)));
% errfro = load(strcat(path, sprintf('errfro_%s.mat',tag)));
% errfro_stream = load(strcat(path, sprintf('errfro_%s-stream.mat',tag)));
% err2_lucp = load(strcat(path, sprintf('err2-base-lucp_%s.mat',tag)));
% errfro_lucp = load(strcat(path, sprintf('errfro-base-lucp_%s.mat',tag)));
%% spectrum
sigma = target.sigma;
figure()
plot(sigma,'k.-')
xlabel('$k$','interpreter','latex')
ylabel('$\sigma_k$','interpreter','latex')
set(gca,'fontsize', 14);
nfro = sqrt(cumsum(sigma.^2,'reverse'));

%% Plot
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
bases = {'full', 'svd', 'sketch'};

% cpqr, fro
err = errfro_cpqr;
optimal = nfro;
subplot(2,2,1)
semilogy(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for idx = 1:length(bases)
    plot(ranks, err.(bases{idx})./optimal(1), strcat(markers{idx},'-'), 'LineWidth', 1.5)
end
hold off
xlim([ranks(1) ranks(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    bases{:},...
    'interpreter','latex')
title(sprintf('CPQR Frobenius norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)

% cpqr, 2
err = err2_cpqr;
optimal = sigma;
subplot(2,2,2)
semilogy(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for idx = 1:length(bases)
    plot(ranks, err.(bases{idx})./optimal(1), strcat(markers{idx},'-'), 'LineWidth', 1.5)
end
hold off
xlim([ranks(1) ranks(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_1$',...
    bases{:},...
    'interpreter','latex')
title(sprintf('CPQR spectral norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)

% lucp, fro
err = errfro_lucp;
optimal = nfro;
subplot(2,2,3)
semilogy(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for idx = 1:length(bases)
    plot(ranks, err.(bases{idx})./optimal(1), strcat(markers{idx},'-'), 'LineWidth', 1.5)
end
hold off
xlim([ranks(1) ranks(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    bases{:},...
    'interpreter','latex')
title(sprintf('LUCP Frobenius norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)

% cpqr, 2
err = err2_lucp;
optimal = sigma;
subplot(2,2,4)
semilogy(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for idx = 1:length(bases)
    plot(ranks, err.(bases{idx})./optimal(1), strcat(markers{idx},'-'), 'LineWidth', 1.5)
end
hold off
xlim([ranks(1) ranks(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_1$',...
    bases{:},...
    'interpreter','latex')
title(sprintf('LUCP spectral norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)
