% DEIM v.s. CPQR v.s. LUPP v.s. RLUPP
clear; close;

nvec = [2000,4000,8000,16000,32000,64000,128000,256000,512000,1024000];
% nvec = [1000, 2000,4000,8000,16000];
% nvec = (2:2:20)*1000;
% kvec = [10,100,500,1000];
kvec = [10,20,30,40];
tag = 'n1024e3_k40';
% tag = 'n1024e3_k1e3';

% time = test_CUR_time(nvec, kvec, tag);
gpu = 1;
time = test_CUR_time(nvec, kvec, tag, gpu);
%%
time = load(sprintf('time_%s.mat',tag));
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
algos = {'LUPP','CPQR','DEIM'};

nplots = length(time.kvec);
figure()
for ik = 1:nplots
    subplot(1,nplots,ik)
    semilogy(time.nvec, time.(algos{1})(:,ik)'./time.nvec, ...
        strcat(markers{1},'-'), 'LineWidth', 1.5) 
    hold on
    for aux = 2:length(algos)
        plot(time.nvec, time.(algos{aux})(:,ik)'./time.nvec, ...
            strcat(markers{aux},'-'), 'LineWidth', 1.5)
    end
    hold off
    xlim([time.nvec(1) time.nvec(end)])
    xlabel('$n$','interpreter','latex')
    ylabel('time / n','interpreter','latex')
    title(sprintf('$k=%d$',time.kvec(ik)), 'interpreter','latex')
    if ik==1, legend(algos{:}, 'interpreter','latex'); end
    set(gca,'FontSize',12)
end