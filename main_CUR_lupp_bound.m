%% load input data
clear; close;
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\dataset\\';
path = 'C:\\Users\\yijundong\\Documents\\MATLAB\\OdenUT\\RandNLA\\CUR\\';
% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'yaleface-64x64';
% tag = 'mnist-train';
% tag = 'GSE10072';
tag = 'ACTIVSg2000';
% tag = 'p2p-Gnutella09';
% tag = 'weightedlaplacian-n1e3-m4e3';
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
k = load(sprintf('rank_%s.mat',tag)); k = k.k;
% test_CUR_rank(k, target, tag, [])
% scp yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/CUR/*_ACTIVSg2000.mat

% tag = strcat(tag,'-srcur');
% algos = {'SRCUR','CPQR','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
% test_CUR_rank(k, target, tag, algos);
%%
l = k + 5;
eta2 = zeros(size(l));
etaf = zeros(size(l));
repeat = 10;

for i =1:length(l)
    eta2sum = 0;
    etafsum = 0;
    for t = 1:repeat
        S = embed(size(A,2), l(i), 'gauss');
        Y = S(A')';
        [L,~,~] = lu(Y,'vec');
        L2 = L(l(i)+1:end,:);
        L1 = L(1:l(i),:);
        aux = (L1'\L2')';
        eta2sum = eta2sum + norm(aux);
        etafsum = etafsum + norm(aux,'fro');
    end
    eta2(i) = eta2sum / repeat;
    etaf(i) = etafsum / repeat;
    fprintf('%d / %d \n', i, length(l))
end
%%
plot(l, eta2)
hold on
plot(l, etaf)
hold off
legend('2','fro')