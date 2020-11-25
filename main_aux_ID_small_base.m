addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\lapack')
%% {DetCPQR, LUCP} for {'full', 'svd', 'sketch'} for column ID
clear; close;
gauss_seed = {'gauss', 500, 500, 500, @(x) log(x)}; % {'gauss', m, n, k, f}
% gauss_seed = {'gauss', 500, 500, 500, @(x) x.^2}; % {'gauss', m, n, k, f}
% snn_seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
% lap_seed = {'laplacian', floor(1e5), floor(5e5)}; % {'laplacian',n,m}
target = TargetMatGenerator(gauss_seed{:});
tag = 'small-gaussian-m500-n500-r500-logdecay';
% tag = 'small-gaussian-m500-n500-r500-polydecay';
% tag = 'small-laplacian-n500-m1000';
plot(target.sigma, 'k.-')

A = target.A;
fprintf('%s \n', tag)
ranks = 20:20:400;
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
% full
base = 'full';
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    
    [j,V] = DetColumnID(A,k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C,V);
    err2_cpqr.(base)(idx) = norm(E);
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% svd
base = 'svd';
W = target.V; % (n,r), r=min(m,n)
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    
    [j,~] = DetColumnID(W(:,1:k)',k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C);
    err2_cpqr.(base)(idx) = norm(E);
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% sketch
base = 'sketch';
os = 10;
m = size(A,1);
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    S = embed(m, k+os);
    X = S(A); % (k+os, n)
    [j,~] = DetColumnID(X,k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C);
    err2_cpqr.(base)(idx) = norm(E);
    errfro_cpqr.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end
%% LUCP
% full
base = 'full';
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    
    [~,j] = CUR_LUCP(A,k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C);
    err2_lucp.(base)(idx) = norm(E);
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% svd
base = 'svd';
W = target.V; % (n,r), r=min(m,n)
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    
    [j,~] = CUR_LUCP(W(:,1:k),k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C);
    err2_lucp.(base)(idx) = norm(E);
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

% sketch
base = 'sketch';
os = 10;
m = size(A,1);
k = ranks(end);
tic;

t = toc;
fprintf('%s for rank %d takes: %.2e sec \n', base, k, t)

for idx = 1:length(ranks)
    k = ranks(idx);
    S = embed(m, k+os);
    X = S(A); % (k+os, n)
    
    [~,j] = CUR_LUCP(X,k);
    
    C = A(:,j(1:k));
    E = ID_Error(A,C);
    err2_lucp.(base)(idx) = norm(E);
    errfro_lucp.(base)(idx) = norm(E,'fro');
    fprintf('%s: %d / %d done \n', base, idx, length(ranks))
end

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
plot(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
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
plot(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
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
plot(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
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
plot(ranks, optimal(ranks+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
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
