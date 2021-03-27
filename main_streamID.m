clear; close;
path_target = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/dataset/';
path = '/Users/ydong/Documents/MATLAB/OdenUT/RandNLA/CUR/';
tag = 'large';
% tag = 'snn-1e3-1e3_a2b1_k100_r1e3_s1e-3';
% tag = 'yaleface-64x64';
% tag = 'mnist-train';
target = load(strcat(path_target, sprintf('target_%s.mat',tag)));
A = target.A;
k = load(sprintf('rank_%s.mat',tag)); k = k.k;
sig = target.sigma;
if target.r == min(size(A)), fprintf('Full rank \n'); end
%%
algos = {'cpqr','lupp','cpqrt','luppt'};
err2 = struct('cpqr', zeros(size(k)), ...
              'lupp', zeros(size(k)), ...
              'cpqrt', zeros(size(k)), ...
              'luppt', zeros(size(k)));
errf = struct('cpqr', zeros(size(k)), ...
              'lupp', zeros(size(k)), ...
              'cpqrt', zeros(size(k)), ...
              'luppt', zeros(size(k)));
            
for i = 1:2
    alg = algos{i};
    for j = 1:length(k)
        [J,V,eta] = CID_stream(A, k(j), alg);
        fprintf('%s : %2d / %2d : %4.2f \n', alg, j, length(k), eta)
        E = A - A(:,J)*V';
        errf.(alg)(j) = norm(E,'fro');
        if issparse(E)
            err2.(alg)(j) = normest(E);
        else
            err2.(alg)(j) = norm(E);
        end
    end
end

for i = 3:4
    alg = algos{i};
    for j = 1:length(k)
        [J,V,eta] = CID_stream(A', k(j), alg);
        fprintf('%s : %2d / %2d : %4.2f \n', alg, j, length(k), eta)
        E = A' - A(J,:)'*V';
        errf.(alg)(j) = norm(E,'fro');
        if issparse(E)
            err2.(alg)(j) = normest(E);
        else
            err2.(alg)(j) = norm(E);
        end
    end
end


%%
save('err2_svcid_large.mat','-struct','err2')
save('errfro_svcid_large.mat','-struct','errf')
err2 = load(fullfile(path,'err2_svcid_large.mat'));
errf = load(fullfile(path,'errfro_svcid_large.mat'));
%%
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
mkmap = struct('cpqr', strcat('r', markers{1},'-'), ...
               'cpqrt', strcat('r', markers{1},':'), ...
               'lupp', strcat('b', markers{2},'-'), ...
               'luppt', strcat('b', markers{2},':'));
markers = arrayfun(@(i) mkmap.(algos{i}), 1:length(algos), 'UniformOutput',false);
legmap = struct('cpqr','CPQR',...
                'lupp','LUPP',...
                'cpqrt','CPQR (m>n)',...
                'luppt','LUPP (m>n)');
labels = arrayfun(@(i) legmap.(algos{i}), 1:length(algos), 'UniformOutput',false);

sigma = target.sigma;
sfro = sqrt(cumsum(sigma.^2,'reverse'));

% frobenius norm
err = errf;
optimal = sfro;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:2
    semilogy(k, (err.(algos{aux}))./optimal(1), markers{aux}, 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_F/||A||_F$','interpreter','latex')
legend('$\sqrt{\sum_{i=k+1}^r \sigma_i^2}/\sqrt{\sum_{i=1}^r \sigma_i^2}$',...
       labels{:},...
       'interpreter','latex')
set(gca,'FontSize',22)

% spectral norm
err = err2;
optimal = sigma;
figure()     
semilogy(k, optimal(k+1)./optimal(1), 'k.-', 'MarkerSize',20, 'LineWidth', 1.5)
hold on
for aux = 1:2
    semilogy(k, (err.(algos{aux}))./optimal(1), markers{aux}, 'LineWidth', 1.5)
end
hold off
xlim([k(1) k(end)])
xlabel('$k$','interpreter','latex')
ylabel('$||A-CUR||_2/||A||_2$','interpreter','latex')
set(gca,'FontSize',22)


