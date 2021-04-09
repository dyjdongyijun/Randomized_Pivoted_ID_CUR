function test_CUR_large(A, ranks, algos, tag)
%%
%   A: (m,n) matrix
%   k: true rank of decomposition (i.e., denoted as l = target rank + O(1) in the manuscript)

    %%
    if isempty(A)
        % default target: (1e6, 1e6) sparse non-negative matrix of rank r = 1000
        m = 1e6; n = m; r = 1000; s = 2/r;
        r0 = 100; a = 2; b = 1;
        X = sprand(m,r,s); Y = sprand(r,n,s); 
        d = 1./(1:r);
        d(1:r0) = d(1:r0)*a; d(r0+1:r) = d(r0+1:r)*b;
        A = (X.*d)*Y;
        
        tag = 'snn-1e6-1e6-a2b1-k100-r1e3-s2e-3';
    end
    
    if ~exist('algos','var') || isempty(algos)
        algos = {'SRCUR',...
                 'CPQR',...
                 'CPQR2pass',...
                 ...
                 'LUPP',...
                 'LUPP2pass',...
                 ...
                 'RSVDDEIM',...
                 ...
                 'RSVDLS'};
    end
    
    if ~exist('tag','var') || isempty(tag)
        tag = 'aux';
    end

    %% Initialization
    rank_file = sprintf('rank_%s.mat',tag);
    if isfile(rank_file)
        k = load(rank_file); k = k.k;
        warning('Rank file exists. Using existing ranks.')
    else
        k = ranks;
    end
    save(rank_file,'k')
    fprintf('write out: %s \n', rank_file)
    
    err_file = sprintf('errfro_%s.mat', tag);
    if isfile(err_file)
        errfro = load(err_file);
    else
        errfro = struct();
    end
    
    time_file = sprintf('time_%s.mat',tag);
    if isfile(time_file)
        time = load(sprintf('time_%s.mat',tag));
    else
        time = struct();
    end
    
    %% Experiments
    for idx = 1:length(algos)
        algo = algos{idx};
        fprintf('%s \n', algo)
        time.(algo) = zeros(size(k));
        errfro.(algo) = zeros(size(k));
        
        
        for t = 1:length(k)
            tic;
            [i,j] = str2curalgo(algo, target, k(t));
            time.(algo)(t) = toc;
            
            fprintf('k = %d: %.4f\n', k(t), time.(algo)(t))
        end

        
        for t = 1:length(k)
            C = A(:,j(1:k(t)));
            R = A(i(1:k(t)),:);
            E = CUR_Error(A,C,R);
            errfro.(algo)(t) = norm(E, 'fro');
            fprintf('%d / %d\t', t, length(k));
        end
        fprintf('\n')
        
        
        save_time(time, tag);
        save_err(errfro, tag)
    end

end



function save_time(time, tag)
    time_file = sprintf('time_%s.mat',tag);
    save(time_file,'-struct','time')
    fprintf('write out: %s \n', time_file)
end

function save_err(err, tag)
    err_file = sprintf('errfro_%s.mat',tag);
    save(err_file,'-struct','err')
    fprintf('write out: %s \n', err_file)
end

