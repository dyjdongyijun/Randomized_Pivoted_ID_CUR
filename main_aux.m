time_inv = 0;
time_ls = 0;
for n = 1000:500:5000
    aux = randn(n);
    
    tic;
    inv(aux);
    time_inv = time_inv + toc;
    
    tic;
    aux\eye(n);
    time_ls = time_ls + toc;
end
fprintf('inv time = %.f, least sqaure time = %.f \n', time_inv, time_ls)
%%
clear; close;
gauss_seed = {'gauss', 500, 500, 500, @(x)x.^(1)}; % {'gauss', m, n, k, f}
snn_seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
lap_seed = {'laplacian', 1000, 5000}; % {'laplacian',n,m}
target = TargetMatGenerator(snn_seed{:});
A = target.A;
plot(target.sigma, 'k.-')
%%
ks = 10:5:50;
randmat = 'gauss';
err_deim = zeros(size(ks));
err_cpqr = zeros(size(ks));
for i = 1:length(ks)
    k = ks(i);
    
    [p_row, p_col] = CUR_ID(A,k,randmat);
    C = A(:,p_col(1:k)); 
    R = A(p_row(1:k),:);
    A_sk = A(p_row(1:k), p_col(1:k));
    E = A - C*(A_sk\R);
    err_cpqr(i) = norm(E,'fro');
    
    [p_row, p_col] = CUR_DEIM(A,k,randmat);
    C = A(:,p_col(1:k)); 
    R = A(p_row(1:k),:);
    U = CUR_U_eval(A,C,R);
    E = A - C*U*R;
    err_deim(i) = norm(E,'fro');
end
%%
sig_fro = sqrt(cumsum(target.sigma.^2,'reverse'));

plot(ks, sig_fro(ks+1), 'k.-')
hold on
plot(ks, err_deim, 'rx-')
plot(ks, err_cpqr, 'bx-')
hold off
legend('Optimal',...
       'DEIM',...
       'CPQR')
%% embed
clear; close;
% seed = {'gauss', 500, 500, 500, @(x)x.^(1)}; % {'gauss', m, n, k, f}
% seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
seed = {'laplacian', 1000, 1000}; % {'laplacian',n,m}
target = TargetMatGenerator(seed{:});
A = target.A;
plot(target.sigma, 'k.-')

ks = 20:20:400;
err_sparse = zeros(size(ks));
err_gauss = zeros(size(ks));
for i = 1:length(ks)
    fprintf('%d \n',ks(i))
    S = embed(size(A,2),ks(i),'sparse5');
    G = embed(size(A,2),ks(i));
    
    Y = S(A')';
    [Q,~] = qr(Y,0);
    E = A - Q*(Q'*A);
    if issparse(E)
        err_sparse(i) = normest(E);
    else
        err_sparse(i) = norm(E);
    end
    
    Y = G(A')';
    [Q,~] = qr(Y,0);
    E = A - Q*(Q'*A);
    if issparse(E)
        err_gauss(i) = normest(E);
    else
        err_gauss(i) = norm(E);
    end
end
%%
plot(ks, target.sigma(ks+1), 'k.-')
hold on
plot(ks, err_gauss, 'rx-')
plot(ks, err_sparse, 'bx-')
hold off
%% embedding time
ns = [1,2,4,8,16,32,64,128,256].*1000;
ks = [20,50,500,1000];
time_g = zeros(length(ks),length(ns));
time_sp3 = zeros(size(time_g));
time_sp5 = zeros(size(time_g));
m = 1000;

for i = 1:length(ks)
    k = ks(i);
    for j = 1:length(ns)
        n = ns(j);
        
        tic;
        S = embed(n,k);
        S(randn(n,m));
        time_g(i,j) = toc;
        
        tic;
        S = embed(n,k,'sparse3');
        S(randn(n,m));
        time_sp3(i,j) = toc;
        
        tic;
        S = embed(n,k,'sparse5');
        S(randn(n,m));
        time_sp5(i,j) = toc;
        
    end
end
%%
nplot = length(ks);
for t = 1:nplot
subplot(1,nplot,t)
plot(ns, time_g(t,:), 'k.-')
hold on
plot(ns, time_sp3(t,:), 'r.-')
plot(ns, time_sp5(t,:), 'b.-')
xlim([ns(1) ns(end)])
title(sprintf('k=%d',ks(t)))
end
legend('gauss','sparse3','sparse5')