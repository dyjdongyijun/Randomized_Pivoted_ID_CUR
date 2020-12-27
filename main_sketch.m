addpath("/Users/ydong/Documents/GraphBLAS/GraphBLAS")
clear; close;
%%
n = 1000;
ks = 2.^(2:11); dfix = 5000;
ds = 500 * 2.^(0:9); kfix = 50;
embeds = {'gauss','srft','sparse3'};
%%
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';
tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'yaleface-64x64';
% tag = 'mnist-train'; 
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
%%
times = struct('ks', ks);
err2s = struct('ks', ks);
errfs = struct('ks', ks);
for idx = 1:length(embeds)
    times.(embeds{idx}) = zeros(2,size(ks,2));
    err2s.(embeds{idx}) = zeros(2,size(ks,2));
    errfs.(embeds{idx}) = zeros(2,size(ks,2));
end
%%
k = 4;
for eidx = 1:length(embeds)
        emb = embeds{eidx};
        S = embed(size(A,2), k, emb);
        Y = S(A')';
end

repeat = 10;
for kidx = 1:length(ks)
    k = ks(kidx);
    for eidx = 1:length(embeds)
        emb = embeds{eidx};
        for aux = 1:repeat
            S = embed(size(A,2), k, emb);
            tic;
            Y = S(A')';
            times.(emb)(2,kidx) = toc;
            times.(emb)(1,kidx) = times.(emb)(1,kidx) + times.(emb)(2,kidx);
            [Q,~] = qr(Y,0);
            E = A - Q*(Q'*A);
            err2s.(emb)(2,kidx) = normest(E);
            errfs.(emb)(2,kidx) = norm(E,'fro');
            err2s.(emb)(1,kidx) = err2s.(emb)(1,kidx) + err2s.(emb)(2,kidx);
            errfs.(emb)(1,kidx) = errfs.(emb)(1,kidx) + errfs.(emb)(2,kidx);
        end
        err2s.(emb)(1,kidx) = err2s.(emb)(1,kidx)/repeat;
        errfs.(emb)(1,kidx) = errfs.(emb)(1,kidx)/repeat;
    end
    fprintf('%d / %d \n', kidx, length(ks))
end
%%
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
legmap = struct('gauss','Gauss',...
                'srft','SRFT',...
                'sparse3','Sparse sign $\zeta=3$');
labels = arrayfun(@(i) legmap.(embeds{i}), 1:length(embeds), 'UniformOutput',false);
%% time - error
xdata = times;
ydata = err2s;
figure()
e = embeds{1};
loglog(xdata.(e)(1,:), ydata.(e)(1,:), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
hold on 
for eidx = 2:length(embeds)
    e = embeds{eidx};
    plot(xdata.(e)(1,:), ydata.(e)(1,:), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
xlabel('time (sec)', 'interpreter', 'latex')
set(gca,'fontsize',16)
%% err-rank v.s. time - rank
sig2 = target.sigma;
sigf = sqrt(cumsum(sig2.^2,'reverse'));

xdata = times;
ydata = err2s;
figure()
e = embeds{1};
loglog(xdata.(e)(1,:), ydata.(e)(1,:), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
hold on 
for eidx = 2:length(embeds)
    e = embeds{eidx};
    plot(xdata.(e)(1,:), ydata.(e)(1,:), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
xlabel('time (sec)', 'interpreter', 'latex')
set(gca,'fontsize',16)
%% subplot(err2, errf, time)
figure()


    
