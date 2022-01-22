function test_CUR_stream(ranks, target, tag, errmeas)
% Input
%   ranks = array of testing ranks
%   target = [struct] fieldnames = {
%                     'A' = m*n matrix,
%                     'description' = matrix description + size
%                     'r' = numerical rank: sigma(r+1)<eps*sigma(1)<=sigma(r)
%                     'sigma' = [min(m,n)*1] singular values
%                     'U','V' = left/right singular vectors}
%   tag = [string] tag for the output files
%   algos = cell of strings representing CUR algorithms being tested on
%   errmeas = [set] {2, 'fro'}
    if ~exist('errmeas','var') || isempty(errmeas)
        errmeas = {2,'fro'};
    end
    
    if ~exist('tag','var') || isempty(tag)
        tag = 'aux';
    end
    %% Initialization
    rank_file = sprintf('rank_%s.mat',tag);
    if isfile(rank_file)
        k = load(rank_file); k = k.k;
        warning('Rank file exists. Using existing ranks.')
        fprintf('Use existing: %s \n', rank_file)
    else
        k = ranks;
        save(rank_file,'k')
        fprintf('write out: %s \n', rank_file)
    end
    
    errs = cell(1, length(errmeas));
    for aux = 1:length(errmeas)
        if isnumeric(errmeas{aux}) && errmeas{aux}==2
            err_file = sprintf('err2_%s.mat', tag);
            if isfile(err_file)
                errs{aux} = load(err_file);
            else
                errs{aux} = struct();
            end
        else
            err_file = sprintf('err%s_%s.mat', errmeas{aux}, tag);
            if isfile(err_file)
                errs{aux} = load(err_file);
            else
                errs{aux} = struct();
            end
        end   
    end
    
    %% Experiments
    algos = {'CPQRstreamCUR', 'CPQRstream'};
    % Initialization
    for t = 1:length(algos)
        algo = algos{t};
        for aux = 1:length(errmeas)
            errs{aux}.(algo) = zeros(size(ranks));
        end
    end
    
    % loop through all ranks
    for t = 1:length(k)
        % run CPQRstreamCUR
        randmat = 'gauss';
        [i,j,U] = CUR_ID_streaming(target.A, k(t), randmat);
        
        % error: approx U
        algo = algos{1};
        A = target.A;
        C = A(:,j(1:k(t)));
        R = A(i(1:k(t)),:);
        for aux = 1:length(errmeas)
            E = A-C*U*R;
            if isnumeric(errmeas{aux}) && errmeas{aux}==2 && issparse(E)
                if max(size(E)) < 5000
                    errs{aux}.(algo)(t) = norm(full(E));
                else
                    errs{aux}.(algo)(t) = normest(E);
                end
            else
                errs{aux}.(algo)(t) = norm(E, errmeas{aux});
            end
        end
        
        fprintf('Approx k = %d: %.4f \n', errs{1}.(algo)(t))
        
        % error: optimal U
        algo = algos{2};
        C = A(:,j(1:k(t)));
        R = A(i(1:k(t)),:);
        E = CUR_Error(A,C,R);
        for aux = 1:length(errmeas)
            if isnumeric(errmeas{aux}) && errmeas{aux}==2 && issparse(E)
                if max(size(E)) < 5000
                    errs{aux}.(algo)(t) = norm(full(E));
                else
                    errs{aux}.(algo)(t) = normest(E);
                end
            else
                errs{aux}.(algo)(t) = norm(E, errmeas{aux});
            end
        end
        
        fprintf('Optimal k = %d: %.4f \n', errs{1}.(algo)(t))
    end

    save_errs(errs, errmeas, tag)
end

%%  % % % % % % % % % % % saving files % % % % % % % % %

function save_errs(errs, errmeas, tag)
    for aux = 1:length(errmeas)
        err = errs{aux};
        if errmeas{aux}==2
            err_file = sprintf('err2_%s.mat', tag);
        else
            err_file = sprintf('err%s_%s.mat', errmeas{aux}, tag);
        end
        save(err_file,'-struct','err')
        fprintf('write out: %s \n', err_file)
    end
end          