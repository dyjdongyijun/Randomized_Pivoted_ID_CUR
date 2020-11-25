%% yijun@natt.oden.utexas.edu
% scp <file> yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/dataset/
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_result_cache = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/result_cache/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
% addpath('/h2/yijun/Documents/MATLAB/lapack')
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\lapack')
addpath('C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\Eigenface')
%% target generation
clear; close;
path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
% path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
% path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
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
path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
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
% test_CUR_rank(k, target, tag, algos);

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
algos = {'SVDDEIM',...
         ...'DetLUPP',...
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
% %% load target
% clear; close;
% target = load('target_GSE10072.mat');
% % A = target.A;
% %% 
% k = 10:10:100;
% err2 = struct('SVDDEIM',zeros(size(k)),...
%               'RSVDDEIM',zeros(size(k)),...
%               'CSSVDDEIM',zeros(size(k)),...
%               'LUCP',zeros(size(k)),...
%               'ACA',zeros(size(k)),...
%               'CSLUCP',zeros(size(k)),...
%               'DetCPQR',zeros(size(k)),...
%               'CPQR',zeros(size(k)),...
%               'CSCPQR',zeros(size(k)),...
%               'CSCPQRstream',zeros(size(k)),...
%               'SVDLS',zeros(size(k)),...
%               'RSVDLS',zeros(size(k)),...
%               'CSSVDLS',zeros(size(k)));
% errfro = struct('SVDDEIM',zeros(size(k)),...
%               'RSVDDEIM',zeros(size(k)),...
%               'CSSVDDEIM',zeros(size(k)),...
%               'LUCP',zeros(size(k)),...
%               'ACA',zeros(size(k)),...
%               'CSLUCP',zeros(size(k)),...
%               'DetCPQR',zeros(size(k)),...
%               'CPQR',zeros(size(k)),...
%               'CSCPQR',zeros(size(k)),...
%               'CSCPQRstream',zeros(size(k)),...
%               'SVDLS',zeros(size(k)),...
%               'RSVDLS',zeros(size(k)),...
%               'CSSVDLS',zeros(size(k)));
% time = struct('SVDDEIM',zeros(size(k)),...
%               'RSVDDEIM',zeros(size(k)),...
%               'CSSVDDEIM',zeros(size(k)),...
%               'LUCP',zeros(size(k)),...
%               'ACA',zeros(size(k)),...
%               'CSLUCP',zeros(size(k)),...
%               'DetCPQR',zeros(size(k)),...
%               'CPQR',zeros(size(k)),...
%               'CSCPQR',zeros(size(k)),...
%               'CSCPQRstream',zeros(size(k)),...
%               'SVDLS',zeros(size(k)),...
%               'RSVDLS',zeros(size(k)),...
%               'CSSVDLS',zeros(size(k))); 
%% Experiments: deterministic
% 
% % SVDDEIM
% fprintf('%s\n', 'SVD-DEIM')
% rsvd = 0;
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_DEIM(A, k(t), rsvd, target.U, target.V);
%     time.SVDDEIM(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.SVDDEIM(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.SVDDEIM(t) = normest(E); 
%     errfro.SVDDEIM(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% % LUCP
% fprintf('%s \n', 'LUCP')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_LUCP(A, k(t));
%     time.LUCP(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.LUCP(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.LUCP(t) = normest(E); 
%     errfro.LUCP(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% % DetCPQR
% rsvd = 0;
% fprintf('%s \n', 'DetCPQR')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_ID(A, k(t), rsvd);
%     time.DetCPQR(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.DetCPQR(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.DetCPQR(t) = normest(E); 
%     errfro.DetCPQR(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% % SVDLS
% rsvd = 'full';
% fprintf('%s \n', 'SVD-LS')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_LeverageScore(A, k(t), rsvd, target.U, target.V);
%     time.SVDLS(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.SVDLS(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.SVDLS(t) = normest(E); 
%     errfro.SVDLS(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% %% Experiments: random
% 
% % RSVDDEIM
% fprintf('%s\n', 'RSVD-DEIM')
% rsvd = 'gauss';
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_DEIM(A, k(t), rsvd);
%     time.RSVDDEIM(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.RSVDDEIM(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.RSVDDEIM(t) = normest(E); 
%     errfro.RSVDDEIM(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % ACA
% fprintf('%s \n', 'ACA')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_LocalSketchACA(A,k(t));
%     time.ACA(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.ACA(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.ACA(t) = normest(E); 
%     errfro.ACA(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % CSLUCP
% fprintf('%s \n', 'CS-LUCP')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_CS_LUCP(A, k(t));
%     time.CSLUCP(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CSLUCP(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CSLUCP(t) = normest(E); 
%     errfro.CSLUCP(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % CPQR
% randmat = 'gauss';
% fprintf('%s \n', 'CPQR')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_ID(A, k(t), randmat);
%     time.CPQR(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CPQR(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CPQR(t) = normest(E); 
%     errfro.CPQR(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % RSVDLS
% rsvd = 'gauss';
% fprintf('%s \n', 'RSVD-LS')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_LeverageScore(A, k(t), rsvd);
%     time.RSVDLS(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.RSVDLS(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.RSVDLS(t) = normest(E,1e-2); 
%     errfro.RSVDLS(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % CSSVDDEIM
% fprintf('%s\n', 'CS-SVD-DEIM')
% randmat = 'sparse3';
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_DEIM(A, k(t), randmat);
%     time.CSSVDDEIM(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CSSVDDEIM(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CSSVDDEIM(t) = normest(E); 
%     errfro.CSSVDDEIM(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% % CSCPQR
% randmat = 'sparse3';
% fprintf('%s \n', 'CS-CPQR')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_ID(A, k(t), randmat);
%     time.CSCPQR(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CSCPQR(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CSCPQR(t) = normest(E); 
%     errfro.CSCPQR(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% % CSSVDLS
% randmat = 'sparse3';
% fprintf('%s \n', 'CS-SVD-LS')
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_LeverageScore(A, k(t), randmat);
%     time.CSSVDLS(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CSSVDLS(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CSSVDLS(t) = normest(E,1e-2); 
%     errfro.CSSVDLS(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
% 
% 
% % CPQR-streaming
% fprintf('%s \n', 'CS-CPQR-stream')
% randmat = 'sparse3';
% for t = 1:length(k)
%     tic;
%     [i,j] = CUR_ID_streaming(A, k(t), randmat);
%     time.CSCPQRstream(t) = toc;
%     fprintf('k = %d: %.4f\n', k(t), time.CSCPQRstream(t))
% end
% 
% for t = 1:length(k)
%     C = A(:,j(1:k(t)));
%     R = A(i(1:k(t)),:);
%     E = CUR_Error(A,C,R);
%     err2.CSCPQRstream(t) = normest(E); 
%     errfro.CSCPQRstream(t) = norm(E,'fro');
%     fprintf('%d / %d\t',t,length(k));
% end
% fprintf('\n')
%% save data
% tag = [];
% save(sprintf('rank_%s.mat',tag),'k')
% save(sprintf('time_%s.mat',tag),'-struct','time')
% save(sprintf('err2_%s.mat',tag),'-struct','err2')
% save(sprintf('errfro_%s.mat',tag),'-struct','errfro')

%% Appendix: spectral + fro norm error plots
% Plot: Deterministic: LUCP/DEIM/CPQR/LS
% sigma = target.sigma;
% sfro = sqrt(cumsum(sigma.^2,'reverse'));

% subplot(1,2,1)
% semilogy(k, sigma(k+1)./sigma(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
% hold on
% semilogy(k, (err2.LUCP)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.SVDDEIM)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.DetCPQR)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.SVDLS)./sigma(1), 's-', 'LineWidth', 1.5)
% hold off
% xlabel('$k$','interpreter','latex')
% ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
% legend('$\sigma_{k+1}/\sigma_1$',...
%     'Deterministic LUCP',...
%     'SVD-DEIM',...
%     'Deterministic CPQR',...
%     'SVD-Leverage score',...
%     'interpreter','latex')
% title(strcat(target.description, '  spectral norm'),...
%     'interpreter','latex')
% set(gca,'FontSize',12)
% 
% subplot(1,2,2)
% semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
% hold on
% semilogy(k, (errfro.LUCP)./sfro(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (errfro.SVDDEIM)./sfro(1), 'o-', 'LineWidth', 1.5)
% semilogy(k, (errfro.DetCPQR)./sfro(1), '^-', 'LineWidth', 1.5)
% semilogy(k, (errfro.SVDLS)./sfro(1), 'd-', 'LineWidth', 1.5)
% hold off
% xlabel('$k$','interpreter','latex')
% ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
% legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
%     'Det-LUCP',...
%     'SVD-DEIM',...
%     'Det-CPQR',...
%     'SVD-LS',...
%     'interpreter','latex')
% title(strcat(target.description, '  Frobenius norm'),...
%     'interpreter','latex')
% set(gca,'FontSize',12)
% 
% 
% Plot: with randomization: DEIM/LUCP(det)/ACA/CSLUCP/CPQR/RSVDLS
% subplot(1,2,1)
% semilogy(k, sigma(k+1)./sigma(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
% hold on
% % semilogy(k, (err2.RSVDDEIM)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.CSSVDDEIM)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.LUCP)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.ACA)./sigma(1), 'd-', 'LineWidth', 1.5)
% semilogy(k, (err2.CSLUCP)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.CPQR)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.CSCPQRstream)./sigma(1), 's-', 'LineWidth', 1.5)
% % semilogy(k, (err2.RSVDLS)./sigma(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (err2.CSSVDLS)./sigma(1), 's-', 'LineWidth', 1.5)
% hold off
% xlabel('$k$','interpreter','latex')
% ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
% legend('$\sigma_{k+1}/\sigma_1$',...
%     ...'RSVD-DEIM',...
%     'CS-SVD-DEIM',...
%     'Deterministic LUCP',...
%     'Local Gaussian Sketch-ACA',...
%     'Count sketch-LUCP',...
%     'Randomized CPQR',...
%     'Streaming CPQR',...
%     ...'RSVD-LS',...
%     'CS-SVD-LS',...
%     'interpreter','latex')
% title(strcat(target.description, '  spectral norm'),...
%     'interpreter','latex')
% set(gca,'FontSize',12)
% 
% subplot(1,2,2)
% semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
% hold on
% semilogy(k, (errfro.RSVDDEIM)./sfro(1), 'o-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CSSVDDEIM)./sfro(1), 'o-', 'LineWidth', 1.5)
% semilogy(k, (errfro.LUCP)./sfro(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (errfro.ACA)./sfro(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CSLUCP)./sfro(1), 's-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CPQR)./sfro(1), '^-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CSCPQR)./sfro(1), '^-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CSCPQRstream)./sfro(1), '^-', 'LineWidth', 1.5)
% semilogy(k, (errfro.RSVDLS)./sfro(1), 'd-', 'LineWidth', 1.5)
% semilogy(k, (errfro.CSSVDLS)./sfro(1), 'd-', 'LineWidth', 1.5)
% hold off
% xlabel('$k$','interpreter','latex')
% ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
% legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
%     'RSVD-DEIM',...
%     'CSSVD-DEIM',...
%     'Det-LUCP',...
%     'LGS-ACA',...
%     'CS-LUCP',...
%     'Rand-CPQR',...
%     'CS-CPQR',...
%     'Stream-CPQR',...
%     'RSVD-LS',...
%     'CSSVD-LS',...
%     'interpreter','latex')
% title(strcat(target.description, '  Frobenius norm'),...
%     'interpreter','latex')
% set(gca,'FontSize',12)
