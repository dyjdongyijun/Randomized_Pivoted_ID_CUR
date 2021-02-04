%%
tsrft = [4.182708e6, 6.431875e6, 6.357542e6, 1.273925e7, 3.3797542e7, 5.9305792e7, 1.13075125e8, 2.05027542e8, 4.61100833e8, 1.050737459e9];
trandn = [2.2845792e7, 5.6160125e7, 1.23514583e8, 2.74439875e8, 4.86020458e8, 9.06236917e8, 1.822762625e9, 1.1397916666e10, 1.8689564916e10, 6.2403764833e10];
tsub = [1.1455333e7, 1.2928458e7, 1.399475e7, 4.664e6, 6.042459e6, 1.5268459e7, 1.0008875e7, 1.7145666e7, 2.097667e7, 2.8249959e7];
ms = 2.^(10:19);
l = 1024;
n = 100;

%%
markers = {'o','s','d','^','v','>','<','p','h','+','*','x'};
%%
loglog(ms, trandn, strcat(markers{1},'-'), 'LineWidth', 1.5)
hold on
plot(ms, tsrft, strcat(markers{2},'-'), 'LineWidth', 1.5)
plot(ms,tsub, strcat(markers{3},'-'), 'LineWidth', 1.5)
hold off
xlim([ms(1), ms(end)])
xlabel('Ambience dimension', 'interpreter', 'latex')
ylabel('Time', 'interpreter', 'latex')
title(sprintf('$l = %d$, $n = %d$', l, n), 'interpreter', 'latex')
legend('Gaussian', 'SRFT', 'Sparse Sign $\zeta=8$', 'interpreter','latex')
set(gca,'FontSize',16)

%%