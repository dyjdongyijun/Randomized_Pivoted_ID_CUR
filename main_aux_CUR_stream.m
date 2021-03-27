%% test
clear; close;
%%
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
% tag = 'ACTIVSg2000';
% tag = 'GSE10072';
% tag = 'large';
% tag = 'p2p-Gnutella09';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'weightedlaplacian-n1e3-m4n';
% tag = 'yaleface-64x64';
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
aug = strcat(tag,'-stream');
k = load(sprintf('rank_%s.mat',aug)); k = k.k;
test_CUR_stream(k, target, aug)
%%
algos = {'CPQRstream',...
        'CPQRstreamCUR'};
% path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
% path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
% tag = 'ACTIVSg2000';
% tag = 'GSE10072';
% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'weightedlaplacian-n1e3-m4n';
% tag = 'yaleface-64x64';
% target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
aug = strcat(tag,'-stream');
k = load(sprintf('rank_%s.mat',aug)); k = k.k;
err2 = load(sprintf('err2_%s.mat',aug));
errfro = load(sprintf('errfro_%s.mat',aug));
%% % % % % % % % % % % Plots % % % % % % % % % %
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
legmap = struct('DetCPQR','Det-CPQR',...
                'CPQR','Rand-CPQR',...
                'CSCPQR','CS-CPQR',...
                'CPQRstream','Stream-CPQR',...
                'CPQRstreamCUR','Stream-CPQR-ApproxU',...
                ...
                'DetLUPP','Det-LUPP',...
                'LUPP','Rand-LUPP',...
                'CSLUPP','CS-LUPP',...
                'LUPPstream','Stream-LUPP',...
                ...
                'SVDDEIM','SVD-DEIM',...
                'RSVDDEIM','RSVD-DEIM',...
                'CSSVDDEIM','CSSVD-DEIM',...
                'CSSVDDEIMstream','Stream-SVD-DEIM',...
                ...
                'LUCP','Det-LUCP',...
                'ACA','DenseSketch-ACA',...(stream)
                'CSLUCP','SparseSketch-LUCP',...(stream)
                ...
                'SVDLS','SVD-LS',...
                'RSVDLS','RSVD-LS',...
                'CSSVDLS','CSSVD-LS',...
                'CSSVDLSstream','Stream-SVD-LS');
%% Plot
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));
% tag = '$10^3 \times 10^3$ SNN($a=2,b=1,k=10^2$,rank$=10^3$)';
% tag = 'Random weighted Laplacian $n=1000, m=4n$';
% tag = 'Dense Gaussian matrix $1000 \times 1000$';
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
title(sprintf('\\texttt{%s} \n Optimal v.s. Approximated $U$',tag),...
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
title(sprintf('\\texttt{%s} \n Optimal v.s. Approximated $U$',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

%% temp
algo = 'CSCPQRstreamCUR';
time.(algo) = zeros(size(k));
errfro.(algo) = zeros(size(k));
a = 8; b = 8;
for t = 1:length(k)
    tic;
    [i,j,U] = CUR_ID_streaming(A, k(t), 'sparse3', a*k(t), b*k(t));
    time.(algo)(t) = toc;
    A = target.A;
    C = A(:,j(1:k(t)));
    R = A(i(1:k(t)),:);
    errfro.(algo)(t) = norm(A-C*U*R,'fro');
    fprintf('k = %d: %.4f\n', k(t), time.(algo)(t))
end
%% rename: CSCPQRstreamCUR -> CSCPQRstreamCUR(const)
oldalgo = sprintf('CSCPQRstreamCUR%dx%d', 8, 8);
newalgo = 'CSCPQRstreamCUR';
time.(newalgo) = time.(oldalgo);
time = rmfield(time, oldalgo);
errfro.(newalgo) = errfro.(oldalgo);
errfro = rmfield(errfro, oldalgo);
%% save update
time_file = sprintf('time_%s.mat',tag);
err_file = sprintf('errfro_%s.mat',tag);
save(time_file,'-struct','time')
save(err_file,'-struct','errfro')
fprintf('update: %s %s \n', time_file, err_file)
%% Post-processing plot
% algos = {'CSSVDDEIM',...
%         'CSSVDDEIMstream',...
%         'CSCPQR',...
%         'CSCPQRstream',...
%         newalgo,...
%         'CSSVDLS',...
%         'CSSVDLSstream'};
sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));

subplot(1,2,1)
semilogy(k, sfro(k+1)./sfro(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
semilogy(k, (errfro.CSSVDDEIM)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSSVDDEIMstream)./sfro(1), 's-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSCPQR)./sfro(1), '^-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSCPQRstream)./sfro(1), '^-', 'LineWidth', 1.5)
semilogy(k, (errfro.(newalgo))./sfro(1), '^-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSSVDLS)./sfro(1), 'd-', 'LineWidth', 1.5)
semilogy(k, (errfro.CSSVDLSstream)./sfro(1), 'd-', 'LineWidth', 1.5)
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')

ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
    'CSSVDDEIM',...
    'CSSVDDEIMstream',...
    'CSCPQR',...
    'CSCPQRstream',...
    newalgo,...
    'CSSVDLS',...
    'CSSVDLSstream',...
    ...'84',...
    'interpreter','latex')
title(sprintf('%s Frobenius norm error',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)

subplot(1,2,2)
plot(k, time.CSSVDDEIM, 's-', 'LineWidth', 1.5)
hold on
plot(k, time.CSSVDDEIMstream, 's-', 'LineWidth', 1.5)
plot(k, time.CSCPQR, '^-', 'LineWidth', 1.5)
plot(k, time.CSCPQRstream, '^-', 'LineWidth', 1.5)
plot(k, time.(newalgo), '^-', 'LineWidth', 1.5)
plot(k, time.CSSVDLS, 'd-', 'LineWidth', 1.5)
plot(k, time.CSSVDLSstream, 'd-', 'LineWidth', 1.5)
hold off
xlim([k(2) k(end)])
ylim([0,2.1])
xlabel('$k$','interpreter','latex')
ylabel('time / sec','interpreter','latex')
legend('CSSVDDEIM',...
    'CSSVDDEIMstream',...
    'CSCPQR',...
    'CSCPQRstream',...
    newalgo,...
    'CSSVDLS',...
    'CSSVDLSstream',...
    'interpreter','latex')
title(sprintf('%s selected runtime',tag),...
    'interpreter','latex')
set(gca,'FontSize',12)