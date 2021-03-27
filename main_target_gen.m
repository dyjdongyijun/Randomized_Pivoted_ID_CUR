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
