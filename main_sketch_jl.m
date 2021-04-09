%% l = 30
randn = [72500.0, 4.369625e6, 3.44675e6, 7.716083e6, 1.142075e7, 2.9093583e7, 6.802225e7, 1.14356292e8, 2.44990292e8, 1.513870792e9];
srft = [86625.0, 1.912833e6, 3.043167e6, 6.194459e6, 1.2449875e7, 2.4662875e7, 4.9789959e7, 9.9380333e7, 1.97864375e8, 4.16063667e8];
sub = [17083.0, 27542.0, 30833.0, 36833.0, 60583.0, 84583.0, 95958.0, 131000.0, 142042.0, 150125.0];
l = 30; n = 100; zeta = 8;
ms = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288];

markers = {'ro','gs','bd','c^','mv','y>','<','p','h','+','*','x'};

figure()
% subplot(1,3,1)
loglog(ms, randn/1e9, strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on
plot(ms, srft/1e9, strcat(markers{2},'-'), 'LineWidth', 1.5)
plot(ms, sub/1e9, strcat(markers{3},'-'), 'LineWidth', 1.5)
hold off
xlim([ms(1), ms(end)])
xlabel('$m$', 'interpreter', 'latex')
ylabel('Time (s)', 'interpreter', 'latex')
title(sprintf('$l = %d$, $n = %d$', l, n), 'interpreter', 'latex')
legend('Gaussian', 'SRFT', sprintf('Sparse Sign $\\zeta=%d$', zeta), 'interpreter','latex')
set(gca,'FontSize',22)

%% l = 100
randn = [487333.0, 6.69925e6, 1.3429917e7, 3.2006959e7, 5.7267625e7, 1.09494958e8, 2.56643708e8, 4.40077791e8, 8.80873292e8, 1.772797708e9];
srft = [297500.0, 1.837083e6, 3.274958e6, 6.183166e6, 1.2772792e7, 2.6718916e7, 5.0638042e7, 1.10509792e8, 2.05587292e8, 4.55185542e8];
sub = [54583.0, 58458.0, 72542.0, 111125.0, 159291.0, 697083.0, 586250.0, 1.4918375e6, 3.565584e6, 7.499833e6];
l = 100; n = 100; zeta = 8;
ms = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288];

figure()
% subplot(1,3,2)
loglog(ms, randn/1e9, strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on
plot(ms, srft/1e9, strcat(markers{2},'-'), 'LineWidth', 1.5)
plot(ms, sub/1e9, strcat(markers{3},'-'), 'LineWidth', 1.5)
hold off
xlim([ms(1), ms(end)])
xlabel('$m$', 'interpreter', 'latex')
ylabel('Time (s)', 'interpreter', 'latex')
title(sprintf('$l = %d$, $n = %d$', l, n), 'interpreter', 'latex')
% legend('Gaussian', 'SRFT', sprintf('Sparse Sign $\\zeta=%d$', zeta), 'interpreter','latex')
set(gca,'FontSize',22)


%% l = 500
randn = [7.267959e6, 3.2516e7, 5.2414667e7, 1.06574209e8, 2.13913875e8, 4.31803375e8, 8.88838666e8, 1.833119208e9, 4.087912292e9, 2.0600064375e10];
srft = [1.954375e6, 3.238e6, 5.361792e6, 1.0346459e7, 2.8071167e7, 5.571775e7, 1.2181575e8, 2.43518958e8, 4.76072792e8, 1.001200208e9];
sub = [1.655209e6, 1.363959e6, 1.683667e6, 1.759291e6, 1.876166e6, 2.226958e6, 2.91725e6, 3.3116417e6, 1.073708e7, 2.3301209e7];
l = 500; n = 100; zeta = 8;
ms = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288];

%% l = 1000
% randn = [2.5412875e7, 5.383825e7, 9.7008167e7, 2.01421042e8, 4.27609667e8, 7.64271083e8, 1.64737075e9, 3.753848208e9, 1.9865646875e10, 5.0470495625e10];
% srft = [4.172666e6, 1.1363208e7, 1.8334667e7, 1.3242625e7, 3.9405417e7, 6.9874792e7, 1.343605e8, 2.64847292e8, 4.82540458e8, 1.140350209e9];
% sub = [1.626e6, 4.996167e6, 8.015875e6, 2.043542e6, 1.896084e6, 2.885167e6, 4.252375e6, 3.7339166e7, 1.954375e6, 7.700208e6];
% l = 1000; n = 100; zeta = 8;
% ms = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288];

%%
figure()
% subplot(1,3,3)
loglog(ms, randn/1e9, strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on
plot(ms, srft/1e9, strcat(markers{2},'-'), 'LineWidth', 1.5)
plot(ms, sub/1e9, strcat(markers{3},'-'), 'LineWidth', 1.5)
hold off
xlim([ms(1), ms(end)])
xlabel('$m$', 'interpreter', 'latex')
ylabel('Time (s)', 'interpreter', 'latex')
title(sprintf('$l = %d$, $n = %d$', l, n), 'interpreter', 'latex')
% legend('Gaussian', 'SRFT', sprintf('Sparse Sign $\\zeta=%d$', zeta), 'interpreter','latex')
set(gca,'FontSize',22)
