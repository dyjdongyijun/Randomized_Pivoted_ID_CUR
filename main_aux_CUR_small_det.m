% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_result_cache = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/result_cache/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_lapack = 'C:\\Users\\yijundong\\Documents\\MATLAB\\lapack';
addpath(path_lapack)
%% target generation
clear; close;
seed = {'gauss', 500, 500, 500, @(x) log(x)}; % {'gauss', m, n, k, f}
% seed = {'gauss', 500, 500, 500, @(x) x.^2}; % {'gauss', m, n, k, f}
% seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
% seed = {'laplacian', floor(500), floor(1000)}; % {'laplacian',n,m}
target = TargetMatGenerator(seed{:});
% tag = 'small-laplacian-n500-m1000';
tag = 'small-gaussian-m500-n500-r500-logdecay';
% tag = 'small-gaussian-m500-n500-r500-polydecay';
plot(target.sigma, 'k.-')
%% test
ranks = 20:20:400;
algos = {'RSVDDEIM',...
         'LUCP',...
         ...'LUCPslow',...
         'DetCPQR',...
         'DetLUPP',...
         'LUPP',...
         'CSLUPP',...
         'LUPPstream'};

errmeas = {2,'fro'};
test_CUR_rank(ranks, target, tag, algos, errmeas)
%% load results
path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
k = load(strcat(path, sprintf('rank_%s.mat',tag))); k = k.k;
time = load(sprintf('time_%s.mat',tag));
errfro = load(sprintf('errfro_%s.mat',tag));
err2 = load(sprintf('err2_%s.mat',tag));
%% % % % % % % % % % % % % % % Plots % % % % % % % % % % % % % % % % % % % 


%% spectrum
sigma = target.sigma;
figure()
plot(sigma,'k.-')
xlabel('$k$','interpreter','latex')
ylabel('$\sigma_k$','interpreter','latex')
set(gca,'fontsize', 14);
nfro = sqrt(cumsum(sigma.^2,'reverse'));

%% errors + runtime
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
% algos = {'SVDDEIM',...
%          'DetLUPP',...
%          'LUCP',...
%          'LUCPslow',...
%          'DetCPQR'};
     
% frobenius norm
err = errfro;
optimal = nfro;
subplot(1,3,1)
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
    algos{:},...
    'interpreter','latex')
title(sprintf('Frobenius norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)


% spectral norm
err = err2;
optimal = sigma;
subplot(1,3,2)
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
legend('$\sigma_{k+1}/\sigma_1$',...
    algos{:},...
    'interpreter','latex')
title(sprintf('Spectral norm error'),...
    'interpreter','latex')
set(gca,'FontSize',14)



% time plot
subplot(1,3,3)
hold on
for aux = 1:length(algos)
    plot(k, time.(algos{aux}), strcat(markers{aux},'-'), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('time / sec','interpreter','latex')
legend(algos{:}, 'interpreter','latex')
title(sprintf('Runtime'),...
    'interpreter','latex')
set(gca,'FontSize',14)

                   