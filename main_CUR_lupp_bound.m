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
% tag = 'GSE10072';
% tag = 'ACTIVSg2000';
% tag = 'p2p-Gnutella09';
% tag = 'weightedlaplacian-n1e3-m4e3';
% target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
% A = target.A;
% k = load(sprintf('rank_%s.mat',tag)); k = k.k;
% test_CUR_rank(k, target, tag, [])
% scp yijun@natt.oden.utexas.edu:/h2/yijun/Documents/MATLAB/RandNLA/CUR/*_ACTIVSg2000.mat

% tag = strcat(tag,'-srcur');
% algos = {'SRCUR','CPQR','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
% test_CUR_rank(k, target, tag, algos);
%%
tags = struct('large','large',...
              'snn','snn-1e3-1e3_a2b1_k100_r1e3_s1e-3',...
              'yaleface','yaleface-64x64',... 
              'mnist','mnist-train');
des = struct('large','large: $4282 \times 8617$',...
              'snn','SNN: $1000 \times 1000$',...
              'yaleface','YaleFace: $165 \times 4096$',... 
              'mnist','MNIST training: $60000 \times 784$');
k = 10:10:150;
%%
repeat = 10;
eta2 = struct('rank',k);
etaf = struct('rank',k);
labels = fieldnames(tags);

for idx = 1:length(labels)
    label = labels{idx};
    tag = tags.(label);
    target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
    A = target.A;
    eta2.(label) = zeros(2,size(k,2));
    etaf.(label) = zeros(2,size(k,2));
    
    for i =1:length(k)
        eta2sum = 0;
        etafsum = 0;
        g2sum = 0;
        gfsum = 0;
        for t = 1:repeat
            S = embed(size(A,2), k(i), 'gauss');
            Y = S(A')';
            [L,~,~] = lu(Y,'vec');
            L2 = L(k(i)+1:end,:);
            L1 = L(1:k(i),:);
            aux = (L1'\L2')';
            eta2sum = eta2sum + norm(aux);
            etafsum = etafsum + norm(aux,'fro');
            
            G = randn(size(Y));
            [Lg,~,~] = lu(G,'vec');
            Lg1 = Lg(1:k(i),:); % (k,k)
            Lg2 = Lg(k(i)+1:end,:); % (m-k,k)
            aux = (Lg1'\Lg2')';
            g2sum = g2sum + norm(aux);
            gfsum = gfsum + norm(aux,'fro');
        end
        eta2.(label)(1,i) = eta2sum / repeat;
        etaf.(label)(1,i) = etafsum / repeat;
        eta2.(label)(2,i) = g2sum / repeat;
        etaf.(label)(2,i) = gfsum / repeat;
        fprintf('%d / %d \n', i, length(k))
    end
end
%%
figure()
for idx = 1:length(labels)
    label = labels{idx};
    tag = tags.(label);
    e2 = eta2.(label)(1,:);
    ef = etaf.(label)(1,:);
    g2 = eta2.(label)(2,:);
    gf = etaf.(label)(2,:);
    k = eta2.rank;
    
    subplot(2,2,idx)
    plot(k, e2, 'r.-', 'MarkerSize',20, 'LineWidth', 1.5)
    hold on
    plot(k, ef, 'b.-', 'MarkerSize',20, 'LineWidth', 1.5)
    plot(k, g2, 'rs:', 'MarkerSize',5, 'LineWidth', 1.5)
    plot(k, gf, 'bs:', 'MarkerSize',5, 'LineWidth', 1.5)
    hold off
    title(des.(label), 'interpreter','latex')
    xlabel('$l$', 'interpreter','latex')
    ylabel('$||L_2 L_1^{-1}||$', 'interpreter','latex')
    if idx==1
        legend('$||.||_2$',...
               '$||.||_F$',...
               'Gaussian $||.||_2$',...
               'Gaussian $||.||_F$',...
               'interpreter','latex')
    end
    set(gca,'Fontsize',12)
end