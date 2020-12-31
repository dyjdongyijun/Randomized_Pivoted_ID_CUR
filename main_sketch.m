addpath("/Users/ydong/Documents/GraphBLAS/GraphBLAS")
clear; close;
%%
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';
% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
tag = 'yaleface-64x64';
% tag = 'mnist-train'; 
target = load(fullfile(path_target, sprintf('target_%s',tag)));
%%
n = 1000;
ks = 10:15:164;
embeds = {'gauss','srft','sparse3'};
test_sketch_rank(tag, ks, embeds, path, path_target)
%%
out = load(fullfile(pwd, sprintf('%s_%s', 'sketch-rank', tag)));
%%
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
legmap = struct('gauss','Gauss',...
                'srft','SRFT',...
                'sparse3','Sparse sign $\zeta=3$');
embeds = out.embeds;
labels = arrayfun(@(i) legmap.(embeds{i}), 1:length(embeds), 'UniformOutput',false);
%% time - error
xdata = out.times;
ydata = out.errfs;

figure()
e = embeds{1};
loglog(xdata.(e), ydata.(e), strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on 
for eidx = 2:length(embeds)
    e = embeds{eidx};
    plot(xdata.(e), ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
xlabel('time (sec)', 'interpreter', 'latex')
set(gca,'fontsize',16)
%% err-rank v.s. time - rank
sig2 = out.sigma;
sigf = sqrt(cumsum(sig2.^2,'reverse'));

figure()
xdata = out.ks;
ydata = out.err2s;
semilogy(xdata, sig2(xdata+1), strcat('.','-'), 'LineWidth', 1.5)
hold on 
for eidx = 1:length(embeds)
    e = embeds{eidx};
    plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
xlabel('rank', 'interpreter', 'latex')
set(gca,'fontsize',16)

figure()
xdata = out.ks;
ydata = out.errfs;
semilogy(xdata, sigf(xdata+1), strcat('.','-'), 'LineWidth', 1.5)
hold on 
for eidx = 1:length(embeds)
    e = embeds{eidx};
    plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_F$', 'interpreter', 'latex')
xlabel('rank', 'interpreter', 'latex')
set(gca,'fontsize',16)

figure()
xdata = out.ks;
ydata = out.times;
e = embeds{1};
semilogy(xdata, ydata.(e), strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on 
for eidx = 2:length(embeds)
    e = embeds{eidx};
    plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
legend(labels{:}, 'interpreter', 'latex')
ylabel('time (sec)', 'interpreter', 'latex')
xlabel('rank', 'interpreter', 'latex')
set(gca,'fontsize',16)

%%

    
