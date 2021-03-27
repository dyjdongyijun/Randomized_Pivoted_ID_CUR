%%
% path = '/h2/yijun/Documents/MATLAB/RandNLA/CUR/';
% path_target = '/h2/yijun/Documents/MATLAB/RandNLA/dataset/';
% path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
% path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';

% tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
tag = 'yaleface-64x64';
% tag = 'mnist-train'; 
% target = load(fullfile(path_target, sprintf('target_%s',tag)));
%%
% n = 1000;
% ks = 30:30:339; % snn-1e3-1e3_a2b1_k100_r1e3_s1e-3
% ks = 10:15:164; % yaleface-64x64
% ks = 50:50:713; % mnist-train
% ks = 40:40:400; % large
% embeds = {'gauss','srft','sparse3'};
% repeat = 10;
% test_sketch_rank(tag, ks, embeds, path, path_target, repeat);
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
% xdata = out.times;
% ydata = out.errfs;
% 
% figure()
% e = embeds{1};
% loglog(xdata.(e), ydata.(e), strcat(markers{1},'-'), 'LineWidth', 1.5)
% hold on 
% for eidx = 2:length(embeds)
%     e = embeds{eidx};
%     plot(xdata.(e), ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
% end
% hold off
% legend(labels{:}, 'interpreter', 'latex')
% ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
% xlabel('time (sec)', 'interpreter', 'latex')
% set(gca,'fontsize',16)

%% err-rank
sig2 = out.sigma;
sigf = sqrt(cumsum(sig2.^2,'reverse'));
if strcmpi(tag,'yaleface-64x64')
    over = 5;
else
    over = 10;
end
fs = 20;

figure()
xdata = out.ks;
ydata = out.err2s;
semilogy(xdata, sig2(xdata-over), strcat('.','-'), 'LineWidth', 1.5)
hold on 
for eidx = 1:length(embeds)
    e = embeds{eidx};
    plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
end
hold off
% legend('$\sigma_{k+1}$', labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_2$', 'interpreter', 'latex')
xlabel('rank', 'interpreter', 'latex')
xlim([xdata(1), xdata(end)])
set(gca,'fontsize',fs)

figure()
xdata = out.ks;
ydata = out.errfs;
semilogy(xdata, sigf(xdata-over), strcat('.','-'), 'LineWidth', 1.5)
hold on 
for eidx = 1:length(embeds)
    e = embeds{eidx};
    plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
    
end
hold off
legend('$\sqrt{\sum_{i > k} \sigma_i^2}$', labels{:}, 'interpreter', 'latex')
ylabel('$||A - Q Q^T A||_F$', 'interpreter', 'latex')
xlabel('rank', 'interpreter', 'latex')
xlim([xdata(1), xdata(end)])
set(gca,'fontsize',fs)

% figure()
% xdata = out.ks;
% ydata = out.times;
% e = embeds{1};
% semilogy(xdata, ydata.(e), strcat(markers{1},'-'), 'LineWidth', 1.5)
% hold on 
% for eidx = 2:length(embeds)
%     e = embeds{eidx};
%     plot(xdata, ydata.(e), strcat(markers{eidx},'-'), 'LineWidth', 1.5)
% end
% hold off
% legend(labels{:}, 'interpreter', 'latex')
% ylabel('time (sec)', 'interpreter', 'latex')
% xlabel('rank', 'interpreter', 'latex')
% xlim([xdata(1), xdata(end)])
% set(gca,'fontsize',fs)

%%

    
