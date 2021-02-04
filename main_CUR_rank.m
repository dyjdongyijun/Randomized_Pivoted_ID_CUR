%% yijun@natt.oden.utexas.edu
% scp <file> yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/dataset/
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_result_cache = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/result_cache/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
% addpath('/h2/yijun/Documents/MATLAB/lapack')
%% target generation
clear; close;
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
% seed = {'gauss', 500, 500, 500, @(x) log(x)}; % {'gauss', m, n, k, f}
% seed = {'gauss', 1000, 1000, 1000, @(x) x.^2}; % {'gauss', m, n, k, f}
% seed = {'snn', 1000, 1000, 100, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
% seed = {'laplacian', floor(500), floor(400)}; % {'laplacian',n,m}
% seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
% target = TargetMatGenerator(seed{:});
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'weightedlaplacian-n1e3-m4e3';
% tag = 'small-gaussian-m1e3-n1e3-r1e3-logdecay';
% tag = 'small-gaussian-m1e3-n1e3-r1e3-polydecay';
% file = sprintf('%starget_%s.mat',path_target,tag);
% save(file, '-struct','target')
% plot(target.sigma, 'k.-')

% Externel data
% clear; close;
tag = 'ACTIVSg2000';
% tag = 'large';
% tag = 'p2p-Gnutella09';
% tag = 'snn-1e4-5e3_a2b1_k10_r1e3_s1e-3';
% tag = 'weightedlaplacian-n1e4-m5n';

% tag = 'GSE10072';
% tag = 'weightedlaplacian-n1e3-m4n';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'small-gaussian-m1e3-n1e3-r1e3-polydecay';
% tag = 'yaleface-64x64';
file = tag;
data_in = load(strcat(path_target, file));
Ain = data_in.Problem.A;
target = TargetMatGenerator(Ain, 'svd');
target.description = tag;
A = target.A;
% plot(target.sigma, 'k.-')
save(fullfile(path_target,sprintf('target_%s',tag)), '-struct', 'target')
%% Eigenface 64x64
clear; close;
load('Yale_64x64.mat'); % fields = {fea = (n*d(=h*w)) image database, gnd = (n*1) label}
fea = fea./max(1e-12,max(fea,[],2)); % normalize 'fea' row-wise btw [0,1]
fea = fea - mean(fea,1);
tag = 'yaleface-64x64';
target = TargetMatGenerator(fea, 'svd');
A = target.A;
target.description = tag;
plot(target.sigma, 'k.-')
file = sprintf('%starget_%s.mat',path_target,tag);
save(file, '-struct','target')

%%
semilogy(target.sigma, 'k.-')
k = 10:10:160;
file = sprintf('%srank_%s.mat',path,tag);
save(file,'k');

%% save target
% tag = 'GSE10072';
% tag = 'ACTIVSg2000';
% tag = 'p2p-Gnutella09';
% tag = 'large';
% tag = 'snn-1e4-5e3_a2b1_k10_r1e3_s1e-3';
% tag = 'weightedlaplacian-n1e4-m5n';
% tag = 'weightedlaplacian-n1e3-m4e3';
% tag = 'yaleface-64x64';
% file = sprintf('target_%s.mat',tag);
% save(strcat(path_target, file),'-struct','target')
%% Plot spectra
% clear; close;
% target1 = load(strcat(path_target, 'target_GSE10072.mat'));
% target2 = load(strcat(path_target, 'target_ACTIVSg2000.mat'));
% target3 = load(strcat(path_target, 'target_p2p-Gnutella09.mat'));
% target4 = load(strcat(path_target, 'target_large.mat'));
% target5 = load(strcat(path_target, 'target_snn-1e4-5e3_a2b1_k10_r1e3_s1e-3.mat'));
% larget6 = load(strcat(path_target, 'target_weightedlaplacian-n1e4-m5n.mat'));
%%
subplot(3,2,1)
semilogy(target1.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('\texttt{GSE10072} spectrum', 'interpreter', 'latex')
xlim([1 length(target1.sigma)])
set(gca,'FontSize',12)

subplot(3,2,2)
semilogy(target2.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('\texttt{ACTIVSg2000} spectrum', 'interpreter', 'latex')
xlim([1 length(target2.sigma)])
set(gca,'FontSize',12)

subplot(3,2,3)
semilogy(target3.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('\texttt{p2p-Gnutella09} spectrum', 'interpreter', 'latex')
xlim([1 length(target3.sigma)])
set(gca,'FontSize',12)

subplot(3,2,4)
semilogy(target4.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('\texttt{large} spectrum', 'interpreter', 'latex')
xlim([1 length(target4.sigma)])
set(gca,'FontSize',12)

subplot(3,2,5)
semilogy(target5.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('SNN($a=2$, $b=1$, $k=10$, rank$=10^3$) spectrum', 'interpreter', 'latex')
xlim([1 length(target4.sigma)])
set(gca,'FontSize',12)

subplot(3,2,6)
semilogy(target6.sigma, 'k.-', 'MarkerSize',10, 'LineWidth', 1)
title('Random weighted Laplacian spectrum', 'interpreter', 'latex')
xlim([1 length(target4.sigma)])
set(gca,'FontSize',12)
%% load input data
% clear; close;
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';
% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
tag = 'yaleface-64x64';
% tag = 'mnist-train';
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
k = load(sprintf('rank_%s.mat',tag)); k = k.k;
% test_CUR_rank(k, target, tag, [])
% scp yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/CUR/*_ACTIVSg2000.mat

tag = strcat(tag,'-srcur');
% algos = {'SRCUR','CPQR','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
algos = {'CPQRstream', 'LUPPstream', 'RSVDDEIMstream', 'RSVDLSstream'};
test_CUR_rank(k, target, tag, algos);

%% load output data
time = load(sprintf('time_%s.mat',tag));
errfro = load(sprintf('errfro_%s.mat',tag));
err2 = load(sprintf('err2_%s.mat',tag));
%% % % % % % % % % % % Plots % % % % % % % % % %
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
legmap = struct('DetCPQR','Det-CPQR',...
                'CPQR','Rand-CPQR',...
                'CPQR2pass','Rand-CPQR-1piter',...
                'CPQR2passOrtho','Rand-CPQR-1piter-ortho',...
                'CSCPQR','CS-CPQR',...
                'CPQRstream','Stream-CPQR',...
                'CSCPQRstream','Stream-CPQR',...
                ...
                'DetLUPP','Det-LUPP',...
                'LUPP','Rand-LUPP',...
                'LUPP2pass','Rand-LUPP-1piter',...
                'LUPP2passOrtho','Rand-LUPP-1piter-ortho',...
                'CSLUPP','CS-LUPP',...
                'LUPPstream','Stream-LUPP',...
                ...
                'SVDDEIM','SVD-DEIM',...
                'RSVDDEIM','RSVD-DEIM',...
                'CSSVDDEIM','CSSVD-DEIM',...
                'RSVDDEIMstream','Stream-SVD-DEIM',...
                'CSSVDDEIMstream','Stream-SVD-DEIM',...
                ...
                'LUCP','Det-LUCP',...
                'ACA','DenseSketch-ACA',...(stream)
                'CSLUCP','SparseSketch-LUCP',...(stream)
                ...
                'SVDLS','SVD-LS',...
                'RSVDLS','RSVD-LS',...
                'CSSVDLS','CSSVD-LS',...
                'RSVDLSstream','Stream-SVD-LS',...
                'CSSVDLSstream','Stream-SVD-LS');
%% Plot: Deterministic: LUPP/LUCP/DEIM/CPQR/LS
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));
% tag = '$10^3 \times 10^3$ SNN($a=2,b=1,k=10^2$,rank$=10^3$)';
% tag = 'Random weighted Laplacian $n=1000, m=4n$';
% tag = 'Dense Gaussian matrix $1000 \times 1000$';
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);
%%
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
         'LUPP2pass',...
         'LUPP2passOrtho',...
         ...
         'CPQR',...
         'CPQR2pass',...
         'CPQR2passOrtho',...
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
         'LUPP2pass',...
         'LUPP2passOrtho',...
         ...
         'CPQR',...
         'CPQR2pass',...
         'CPQR2passOrtho',...
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
algos = {'CSSVDDEIMstream',...
         ...
         'LUPPstream',...
         ...
         ...'CSLUCP',...
         'ACA',...
         ...
         'CSCPQRstream',...
         ...
         'CSSVDLSstream'};
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

%% Appendix: experiments

