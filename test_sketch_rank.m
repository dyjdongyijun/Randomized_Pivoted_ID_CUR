function test_sketch_rank(tag, ks, embeds, path, path_large, repeat)
    if ~exist('ks','var') || isempty(ks)
        ks = 10:10:100;
    end
    if ~exist('embeds','var') || isempty(embeds)
        embeds = {'gauss','srft','sparse3'};
    end
    if ~exist('path','var') || isempty(path)
        path = pwd;
    end
    if ~exist('path_large','var') || isempty(path_large)
        path_large = pwd;
    end
    if ~exist('repeat','var') || isempty(repeat)
        repeat = 5;
    end
    out = struct();
    out.tag = tag;
    out.ks = ks;
    out.embeds = embeds;
    out.repeat = repeat;
    %% load target
    target = loader(tag, 'target', path_large);
    A = target.A;
    %% initiate
    times = struct();
    err2s = struct();
    errfs = struct();
    for idx = 1:length(embeds)
        times.(embeds{idx}) = zeros(size(ks));
        err2s.(embeds{idx}) = zeros(size(ks));
        errfs.(embeds{idx}) = zeros(size(ks));
    end
    %% experiments
    k = 4;
    for eidx = 1:length(embeds)
        emb = embeds{eidx};
        embed(size(A,2), k, emb);
    end

    for kidx = 1:length(ks)
        k = ks(kidx);
        for eidx = 1:length(embeds)
            emb = embeds{eidx};
            for aux = 1:repeat
                S = embed(size(A,2), k, emb);
                tic;
                Y = S(A')';
                taux = toc;
                times.(emb)(kidx) = times.(emb)(kidx) + taux;
                [Q,~] = qr(Y,0);
                E = A - Q*(Q'*A);
                err2 = normest(E);
                errf = norm(E,'fro');
                err2s.(emb)(kidx) = err2s.(emb)(kidx) + err2;
                errfs.(emb)(kidx) = errfs.(emb)(kidx) + errf;
            end
            err2s.(emb)(kidx) = err2s.(emb)(kidx)/repeat;
            errfs.(emb)(kidx) = errfs.(emb)(kidx)/repeat;
        end
        fprintf('%d / %d \n', kidx, length(ks))
    end

    %% save
    out.sigma = target.sigma;
    out.times = times;
    out.err2s = err2s;
    out.errfs = errfs;
    name = fullfile(path, sprintf('%s_%s.mat', 'sketch-rank', tag));
    if exist(name,'file')==2
        fprintf('%s exists and updated \n', name)
    end
    save(name, '-struct', 'out');
    fprintf('write out: %s \n', name);
end

%% 
function data = loader(tag, prefix, path)
    if ~exist('path','var') || isempty(path)
        path = pwd;
    end
    name = fullfile(path, sprintf('%s_%s.mat', prefix, tag));
    data = load(name);
end