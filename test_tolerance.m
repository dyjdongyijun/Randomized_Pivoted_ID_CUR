gen = @(n,spec) (randn(n,length(spec)).*spec)*randn(length(spec),n)./length(spec);

n = 10000; spec = (1:n).^(-2);
A = gen(n,spec);
%%
X = randn(floor(n/log2(n)),n)*A;
tic;
[~,~,~] = qr(full(X),'vec');
t = toc;
fprintf('CPQR time = %8.4f sec \n',t)

tic;
s = svd(X);
[~,~,~] = lu(full(X'),'vec');
t = toc;
fprintf('SVD time = %8.4f sec \n',t)
%%
ref = norm(A,'fro');
tol = 1e-1;
fprintf('Tolerance = %4.1e ref = %8.4e \n', tol, tol*ref);
tic;
[J1,VT1] = rand_id(A,tol*ref,[],'lupp');
t = toc;
fprintf('lupp accuracy = %8.4e, time = %8.4f \n',norm(A - A(:,J1)*VT1,'fro')/ref, t)
tic;
[J2,VT2] = rand_id(A,tol*ref,[],'cpqr');
t = toc;
fprintf('cpqr accuracy = %8.4e, time = %8.4f \n',norm(A - A(:,J2)*VT2,'fro')/ref, t)
%% some test results
% 
% accuracy = relative error
% 
% Target matrix
% n = 10000; spec = (1:n).^(-2);
% A = gen(n,spec);
% 
% Tolerance = 1.0e-01 ref = 1.0448e-01 
% lupp accuracy = 9.1979e-03, time =   3.1791 
% cpqr accuracy = 3.8194e-02, time =   3.5128 
% 
% Tolerance = 1.0e-02 ref = 1.0448e-02 
% lupp accuracy = 2.0536e-03, time =   5.3097 
% cpqr accuracy = 3.5269e-03, time =   4.9305 
% 
% Tolerance = 1.0e-03 ref = 1.0448e-03 
% lupp accuracy = 3.8183e-04, time =   3.4156 
% cpqr accuracy = 2.2537e-04, time =   4.4521
% 
% Tolerance = 1.0e-04 ref = 1.0448e-04 
% lupp accuracy = 5.2768e-05, time =   7.0687 
% cpqr accuracy = 5.2682e-05, time =   7.3441 
% 
% Tolerance = 1.0e-04 ref = 1.0448e-04 
% lupp accuracy = 8.3291e-05, time =   5.4333 
% cpqr accuracy = 5.2630e-05, time =   7.2112 
% 
% Tolerance = 1.0e-05 ref = 1.0448e-05 
% lupp accuracy = 5.1499e-05, time =   6.2790 
% cpqr accuracy = 4.9838e-05, time =   7.3612 
% 
% Tolerance = 1.0e-08 ref = 1.0448e-08 
% lupp accuracy = 5.1545e-05, time =   6.4829 
% cpqr accuracy = 4.9972e-05, time =   6.8820 