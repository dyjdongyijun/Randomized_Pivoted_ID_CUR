function time = test_CUR_time(nvec, kvec, tag, gpu)
    if ~exist('nvec','var') || isempty(nvec)
        nvec = (2:2:20)*1000;
    end
    if ~exist('kvec','var') || isempty(kvec)
        kvec = [10,100,500,1000];
    end
    if ~exist('tag','var') || isempty(tag)
        tag = 'aux';
    end
    if ~exist('gpu','var')
        gpu = [];
    end
%%
    time = struct(...
              'tag', tag, ...
              'nvec', nvec,...
              'kvec', kvec,...
              'LUPP', zeros(length(nvec),length(kvec)),...
              'CPQR', zeros(length(nvec),length(kvec)),...
              'DEIM', zeros(length(nvec),length(kvec))...
              );

    for ni = 1:length(nvec)
        n = nvec(ni);
        for ki = 1:length(kvec)
            k = kvec(ki);
            Y= randn(n,k);
            if ~isempty(gpu)
                gY = gpuArray(Y);
            end

            % LUPP
            if isempty(gpu)
                tic
                [~,~,i] = lu(Y,'vector');
                time.LUPP(ni,ki) = toc;
            else
                tic
                [~,~,i] = lu(gY,'vector');
                time.LUPP(ni,ki) = toc;
            end
                

            % DEIM
            if isempty(gpu)
                tic
                [Q,R] = qr(Y,0);
                [Uk,~,~] = svd(R,'econ');
                [~,~,i] = lu(Q*Uk,'vector');
                time.DEIM(ni,ki) = toc;
            else
                tic
                [Q,R] = qr(gY,0);
                [Uk,~,~] = svd(R,'econ');
                [~,~,i] = lu(Q*Uk,'vector');
                time.DEIM(ni,ki) = toc;
            end

            % CPQR
            tic
            [~,~,i] = qr(Y','vector');
            time.CPQR(ni,ki) = toc;

            save_time(time, tag)
            fprintf('n=%6d   k=%4d   t_lupp = %8.3f    t_deim = %8.3f    t_cpqr = %8.3f\n',...
                    n,k,time.LUPP(ni,ki),time.DEIM(ni,ki),time.CPQR(ni,ki))
        end
    end
end
%% 
function save_time(time, tag)
    time_file = sprintf('time_%s.mat',tag);
    save(time_file,'-struct','time')
    fprintf('write out: %s \n', time_file)
end