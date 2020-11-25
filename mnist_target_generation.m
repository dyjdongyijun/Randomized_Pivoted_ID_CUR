% addpath('/h2/yijun/Documents/MATLAB/lapack')
clear; close;
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\lapack')
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR')
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR_mnist')
path_mnist = 'C:\\Users\\yijundong\\Documents\\MATLAB\\dataset\\mnist';
path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR_mnist';
path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
%%
path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR_mnist/';
path_mnist = '/h2/yijun/Documents/MATLAB/dataset/mnist';
addpath('/h2/yijun/Documents/MATLAB/lapack/')
addpath('/h2/yijun/Documents/MATLAB/RandNLA/CUR')
addpath('/h2/yijun/Documents/MATLAB/RandNLA/CUR_mnist')
%% load MNIST
fimgtrain = fullfile(path_mnist, 'train-images.idx3-ubyte');
flabtrain = fullfile(path_mnist, 'train-labels.idx1-ubyte');
fimgtest = fullfile(path_mnist, 't10k-images.idx3-ubyte');
flabtest = fullfile(path_mnist, 't10k-labels.idx1-ubyte');

[img_train, lab_train] = mnist_parse(fimgtrain, flabtrain);
[img_test, lab_test] = mnist_parse(fimgtest, flabtest);
num_train = size(lab_train,1);
num_test = size(lab_test,1);
Atrain = double(reshape(img_train,[],num_train)');
size(Atrain)
Atest = double(reshape(img_test,[],num_test)');
size(Atest)
tags = {'mnist-train','mnist-test'};
As = {Atrain, Atest};
nsmp = [num_train, num_test];
%% Target generation
aux = 1;
tag = tags{aux};
A = As{aux};
m = nsmp(aux);
A = A./max(A,[],'all');
target = TargetMatGenerator(A,'svd');
target.description = tag;
A = target.A;
save(fullfile(path_target,sprintf('target_%s',tag)), '-struct', 'target')
%% save target
semilogy(target.sigma, 'k.-')
% file = sprintf('target_%s.mat',path_target,tag);
% save(file, '-struct','target')
%%
k = 400;
rsvd = 0;
[i,j] = CUR_DEIM(target.A, k, rsvd, target.U, target.V);
%%
lab = lab_train(i);
figure();
histogram(lab)

n_smp = 32;
sample_idx = randperm(m, n_smp);
sp1 = 8; sp2 = 8;

figure();
for t = 1:n_smp
    subplot(sp1,sp2,1+2*(t-1))
    imshow(full(reshape(A(sample_idx(t),:), 28, [])))
    subplot(sp1,sp2,2*t)
    imshow(full(reshape(A(sample_idx(t),j), 20, [])))
end
%% scale rank
k = 40:40:400;
file = fullfile(path, sprintf('rank_%s.mat',tag));
save(file,'k');
%%
k = load(fullfile(path, sprintf('rank_%s.mat',tag))); k = k.k;
%%
test_CUR_rank(k, target, tag, [])    
% cpqrs = {'DetCPQR','CPQR','CSCPQR','CPQRstream'};
% lupps = {'DetLUPP','LUPP','CSLUPP','LUPPstream',};
% deims = {'SVDDEIM','RSVDDEIM','CSSVDDEIM','RSVDDEIMstream'};
% lucps = {'LUCP','ACA','CSLUCP'};
% lss = {'SVDLS','RSVDLS','CSSVDLS','RSVDLSstream'};
% test_CUR_rank(k, target, tag, cpqrs)    
% test_CUR_rank(k, target, tag, lupps)
% test_CUR_rank(k, target, tag, deims)
% test_CUR_rank(k, target, tag, lucps)
% test_CUR_rank(k, target, tag, lss)
% scp yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/CUR/*_ACTIVSg2000.mat
%% load output data
time = load(sprintf('time_%s.mat',tag));
errfro = load(sprintf('errfro_%s.mat',tag));
err2 = load(sprintf('err2_%s.mat',tag));
%% % % % % % % % % % % Plots % % % % % % % % % %
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
legmap = struct('DetCPQR','Det-CPQR',...
                'CPQR','Rand-CPQR',...
                'CSCPQR','CS-CPQR',...
                'CSCPQRstream','Stream-CPQR',...
                ...
                'DetLUPP','Det-LUPP',...
                'LUPP','Rand-LUPP',...
                'CSLUPP','CS-LUPP',...
                'LUPPstream','Stream-LUPP',...
                ...
                'SVDDEIM','SVD-DEIM',...
                'RSVDDEIM','RSVD-DEIM',...
                'CSSVDDEIM','CSSVD-DEIM',...
                'RSVDDEIMstream','Stream-SVD-DEIM',...
                ...
                'LUCP','Det-LUCP',...
                'ACA','DenseSketch-ACA',...(stream)
                'CSLUCP','SparseSketch-LUCP',...(stream)
                ...
                'SVDLS','SVD-LS',...
                'RSVDLS','RSVD-LS',...
                'CSSVDLS','CSSVD-LS',...
                'RSVDLSstream','Stream-SVD-LS');
%% Plot: Deterministic: LUPP/LUCP/DEIM/CPQR/LS
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));
% tag = '$10^3 \times 10^3$ SNN($a=2,b=1,k=10^2$,rank$=10^3$)';
% tag = 'Random weighted Laplacian $n=1000, m=4n$';
% tag = 'Dense Gaussian matrix $1000 \times 1000$';
algos = {'SVDDEIM',...
         'DetLUPP',...
         'LUCP',...
         'DetCPQR',...
         'SVDLS'};
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);

% frobenius norm
err = errfro;
optimal = sfro;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Deterministic',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

% spectral norm
err = err2;
optimal = sigma;
figure()
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_{1}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Deterministic',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

%% Plot: randomized, non-stream
algos = {'RSVDDEIM',...
         ...
         'LUPP',...
         ...
         'CPQR',...
         ...
         'RSVDLS'};
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);

% frobenius norm
err = errfro;
optimal = sfro;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Randomized',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

% spectral norm
err = err2;
optimal = sigma;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_{1}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Randomized',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

%% Plot: randomized, stream
algos = {'RSVDDEIMstream',...
         ...
         'LUPPstream',...
         ...
         'CSLUCP',...
         'ACA',...
         ...
         'CSCPQRstream',...
         ...
         'RSVDLSstream',...
         };
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);

% frobenius norm
err = errfro;
optimal = sfro;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Single-view',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

% spectral norm
err = err2;
optimal = sigma;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_{1}$',...
    labels{:},...
    'interpreter','latex')
title(sprintf('\\texttt{%s} Single-view',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

%% time plot: non-stream
algos = {'RSVDDEIM',...
         ...
         'LUPP',...
         ...
         'CPQR',...
         ...
         'RSVDLS'};
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);
     
figure()
semilogy(k, time.(algos{1}), strcat(markers{1},'-'), 'LineWidth', 1.5) 
hold on
for aux = 2:length(algos)
    plot(k, time.(algos{aux}), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(2) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('time / sec','interpreter','latex')
legend(labels{:}, 'interpreter','latex')
title(sprintf('\\texttt{%s} randomized',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)


% time plot: stream
algos = {'RSVDDEIMstream',...
         ...
         'LUPPstream',...
         ...
         'CSLUCP',...
         'ACA',...
         ...
         'CSCPQRstream',...
         ...
         'RSVDLSstream'};
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);

figure()
semilogy(k, time.(algos{1}), strcat(markers{1},'-'), 'LineWidth', 1.5) 
hold on
for aux = 2:length(algos)
    plot(k, time.(algos{aux}), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(2) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('time / sec','interpreter','latex')
legend(labels{:}, 'interpreter','latex')
title(sprintf('\\texttt{%s} single-view',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)
%%


