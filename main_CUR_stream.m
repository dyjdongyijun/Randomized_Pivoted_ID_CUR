%% yijun@natt.oden.utexas.edu
% scp <file> yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/dataset/
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_result_cache = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/result_cache/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
% addpath('/h2/yijun/Documents/MATLAB/lapack')

%% load input data
clear; close;
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';
% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'yaleface-64x64';
% tag = 'mnist-train';
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
k = load(sprintf('rank_%s.mat',tag)); k = k.k;
% test_CUR_rank(k, target, tag, [])
% scp yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/CUR/*_ACTIVSg2000.mat

% tag = strcat(tag,'-srcur');
% algos = {'SRCUR','CPQR','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
% algos = {'CPQRstream', 'LUPPstream', 'RSVDDEIMstream', 'RSVDLSstream'};
% test_CUR_rank(k, target, tag, algos);

%% load output data
time = load(sprintf('time_%s.mat',tag));
errfro = load(sprintf('errfro_%s.mat',tag));
err2 = load(sprintf('err2_%s.mat',tag));
%% % % % % % % % % % % Plots % % % % % % % % % %
algos = {'RSVDDEIMstream','LUPPstream','CPQRstream','RSVDLSstream'};

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
                'CSSVDLSstream','Stream-SVD-LS',...
                ...
                'SRCUR','SRCUR');
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);


mkmap = struct('RSVDDEIMstream', strcat('r', markers{1},'-'), ...
               'LUPPstream', strcat('b', markers{2},'-'), ...
               'CPQRstream', strcat('g', markers{3},'-'), ...
               'RSVDLSstream', strcat('c', markers{4},'-'));
markers = arrayfun(@(i) mkmap.(algos{i}), 1:length(algos), 'UniformOutput',false);


sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));
% tag = '$10^3 \times 10^3$ SNN($a=2,b=1,k=10^2$,rank$=10^3$)';
% tag = 'Random weighted Laplacian $n=1000, m=4n$';
% tag = 'Dense Gaussian matrix $1000 \times 1000$';

% frobenius norm
err = errfro;
optimal = sfro;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), mkmap.(algos{aux}), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
       labels{:},...
       'interpreter','latex')
% title(sprintf('\\texttt{%s} Randomized',tag), 'interpreter','latex')
set(gca,'FontSize',22)

% spectral norm
err = err2;
optimal = sigma;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:length(algos)
    semilogy(k, (err.(algos{aux}))./optimal(1), mkmap.(algos{aux}), 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
% legend('$\sigma_{k+1}/\sigma_{1}$', labels{:}, 'interpreter','latex')
% title(sprintf('\\texttt{%s} Randomized',tag), 'interpreter','latex')
set(gca,'FontSize',22)

% time
figure()
semilogy(k, time.(algos{1}), mkmap.(algos{aux}), 'LineWidth', 1.5) 
hold on
for aux = 1:length(algos)
    plot(k, time.(algos{aux}), mkmap.(algos{aux}), 'LineWidth', 1.5)
end
hold off
xlim([k(2) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('Time (s)','interpreter','latex')
% legend(labels{:}, 'interpreter','latex')
% title(sprintf('\\texttt{%s} single-view',tag), 'interpreter','latex')
set(gca,'FontSize',22)
