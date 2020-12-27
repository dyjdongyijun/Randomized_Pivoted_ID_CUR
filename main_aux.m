mls = [1000,2000,4000,8000,16000];
kls = 20:20:200;
eta = zeros(length(mls), length(kls));
repeat = 10;
linear = @(m) 1:m;
poly = @(m) (1:m).^3;
expo = @(m) exp(1:m);

for im = 1:length(mls)
    m = mls(im);
    for ik = 1:length(kls)
        k = kls(ik);
        aux = 0;
        for t = 1:repeat
            Y = randn(m,k);
            [L,~,~] = lu(Y,'vec');
            L1 = L(1:k,:); % (k,k)
            L2 = L(k+1:end,:); % (m-k,k)
            aux = aux + norm(L1'\L2');
        end
        eta(im,ik) = aux / repeat;
    end
end
%%
figure()
for im = 1:length(mls)
    m = mls(im);
    subplot(2,2,im)
    plot(kls, eta(im,:), '.-', 'MarkerSize',20, 'LineWidth', 1.5)
    title(sprintf('m=%d', m))
    xlabel('k')
    ylabel('\eta')
end
%%