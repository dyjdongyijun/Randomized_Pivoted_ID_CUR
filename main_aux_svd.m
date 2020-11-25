clear; close;
gauss_seed = {'gauss', 500, 500, 500, @(x)x.^{3}}; % {'gauss', m, n, k, f}
snn_seed = {'snn', 1000, 1000, 20, 1000, 2, 1, 1e-3}; % {'snn',m,n,k,r,a,b,s}
lap_seed = {'laplacian', 1000, 5000}; % {'laplacian',n,m}
target = TargetMatGenerator(lap_seed{:});
A = target.A;
plot(target.sigma, 'k.-')
%%
k = 100; randmat = 'sparse3'; l = 2*k; s = 2*l;
[U_stream,d_stream,V_stream] = stream_svd(A, k, randmat, l, s);
[U,d,V] = RSVD(A,k,randmat);
% [U,D,V] = svd(full(A),0); d = diag(D);
%%
close;
semilogy(target.sigma(1:k),'k.-')
hold on
plot(d_stream, 'rx-')
plot(d,'bx-')
% plot(d(1:k),'bo-')
hold off
title(sprintf('$||U_{stream}^T U||_2 = %4.2f$, $||V_{stream}^T V||_2 = %4.2f$',...
              norm(U_stream'*target.U), norm(V_stream'*target.V)),...
      'interpreter','latex')
ylabel('$\sigma_k$', 'interpreter','latex')
xlabel('$k$','interpreter','latex')
legend('SVD',...
       'Streaming SVD',...
       'RSVD',...
       'interpreter','latex')
%%