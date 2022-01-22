function test_CUR_rank(ranks, target, tag, algos, errmeas)
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
    
    if ~exist('algos','var') || isempty(algos)
        algos = {...'SRCUR',...
                'DetCPQR',...
                'CPQR',...
                'CPQR2pass',...
                'CPQR2passOrtho',...
                ...'CSCPQR',...
                'CPQRstream',...
                'CPQRstreamCUR',...
                ...
                ...'DetLUPP',...
                'LUPP',...
                'LUPP2pass',...
                'LUPP2passOrtho',...
                ...'CSLUPP',...
                'LUPPstream',...
                ...
                ...'SVDDEIM',...
                'RSVDDEIM',...
                ...'CSSVDDEIM',...
                'RSVDDEIMstream',...
                ...
                ...'LUCP',...
                'ACA',...(stream)
                ...'CSLUCP',...(stream)
                ...
                ...'SVDLS',...
                'RSVDLS',...
                ...'CSSVDLS',...
                'RSVDLSstream'...
                };
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
    
    errs = cell(1, length(errmeas));
    for aux = 1:length(errmeas)
        if errmeas{aux}==2
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
    
    time_file = sprintf('time_%s.mat',tag);
    if isfile(time_file)
        time = load(sprintf('time_%s.mat',tag));
    else
        time = struct();
    end
    A = target.A;
    %% Experiments
    for idx = 1:length(algos)
        algo = algos{idx};
        fprintf('%s \n', algo)
        time.(algo) = zeros(size(k));
        for aux = 1:length(errmeas)
            errs{aux}.(algo) = zeros(size(k));
        end
        
        for t = 1:length(k)
            if strcmpi(algo, 'CPQRstreamCUR')
                tic; 
                [i,j,U] = str2curalgo(algo, target, k(t));
                time.(algo)(t) = toc;
                
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
            else
                tic;
                [i,j] = str2curalgo(algo, target, k(t));
                time.(algo)(t) = toc;
            end
            fprintf('k = %d: %.4f\n', k(t), time.(algo)(t))
        end

        if ~strcmpi(algo, 'CPQRstreamCUR')
            for t = 1:length(k)
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
                fprintf('%d / %d\t', t, length(k));
            end
            fprintf('\n')
        end
        
        save_time(time, tag);
        save_errs(errs, errmeas, tag)
    end
    
end

%%  % % % % % % % % % % % saving files % % % % % % % % %
function save_time(time, tag)
    time_file = sprintf('time_%s.mat',tag);
    save(time_file,'-struct','time')
    fprintf('write out: %s \n', time_file)
end

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

%% % % % % % % % % % % % test functions for each algo % % % % % % % % % 
function [i,j,U] = str2curalgo(algo, target, rank)
% % % % % SRCUR % % % % % 
    if strcmpi(algo, 'SRCUR')
        [i,j] = srcur(target.A, rank);
% % % % % LUPP % % % % % 
    elseif strcmpi(algo, 'DetLUPP')
        sketch = [];
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'LUPP')
        sketch = 'gauss';
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'LUPP2pass')
        sketch = 'gauss';
        stream = [];
        power = 1;
        ortho = 0;
        [i,j] = CUR_LUPP(target.A, rank, sketch, stream, power, ortho);
    elseif strcmpi(algo, 'LUPP2passOrtho')
        sketch = 'gauss';
        stream = [];
        power = 1;
        ortho = 1;
        [i,j] = CUR_LUPP(target.A, rank, sketch, stream, power, ortho);
    elseif strcmpi(algo, 'CSLUPP')
        sketch = 'sparse3';
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'LUPPstream')
        sketch = 'gauss';
        stream = 1;
        [i,j] = CUR_LUPP(target.A, rank, sketch, stream);
% % % % % CPQR % % % % % 
    elseif strcmpi(algo, 'DetCPQR')
        sketch = [];
        [i,j] = CUR_ID(target.A, rank, sketch);
    elseif strcmpi(algo, 'CPQR2pass')
        sketch = 'gauss';
        piter = 1;
        ortho = 0;
        [i,j] = CUR_ID(target.A, rank, sketch, piter, ortho);
    elseif strcmpi(algo, 'CPQR2pass')
        sketch = 'gauss';
        piter = 1;
        ortho = 1;
        [i,j] = CUR_ID(target.A, rank, sketch, piter, ortho);
    elseif strcmpi(algo, 'CPQR')
        sketch = 'gauss';
        [i,j] = CUR_ID(target.A, rank, sketch);
    elseif strcmpi(algo, 'CSCPQR')
        sketch = 'sparse3';
        [i,j] = CUR_ID(target.A, rank, sketch);
    elseif strcmpi(algo, 'CPQRstream')
        sketch = 'gauss';
        l = min(min(size(target.A)), rank+10);
        s = l;
        [i,j] = CUR_ID_streaming(target.A, rank, sketch, l, s);
    elseif strcmpi(algo, 'CPQRstreamCUR')
        sketch = 'gauss';
        [i,j,U] = CUR_ID_streaming(target.A, rank, sketch);
% % % % % LUCP % % % % % 
    elseif strcmpi(algo,'LUCP')
        [i,j] = CUR_LUCP(target.A, rank);
    elseif strcmpi(algo,'LUCPslow')
        disable_lapack = 1;
        [i,j] = CUR_LUCP(target.A, rank, disable_lapack);
    elseif strcmpi(algo,'ACA')
        [i,j] = CUR_ACA(target.A, rank);
    elseif strcmpi(algo, 'CSLUCP')
        stream = 1;
        [i,j] = CUR_CS_LUCP(target.A, rank, stream);
% % % % % DEIM % % % % % 
    elseif strcmpi(algo,'SVDDEIM')
        rsvd = 0;
        [i,j] = CUR_DEIM(target.A, rank, rsvd, target.U, target.V);
    elseif strcmpi(algo,'RSVDDEIM')
        rsvd = 'gauss';
        [i,j] = CUR_DEIM(target.A, rank, rsvd);
    elseif strcmpi(algo, 'CSSVDDEIM')
        randmat = 'sparse3';
        [i,j] = CUR_DEIM(target.A, rank, randmat);
    elseif strcmpi(algo, 'RSVDDEIMstream')
        randmat = 'gauss';
        stream = 1;
        [i,j] = CUR_DEIM(target.A, rank, randmat, stream);
% % % % % LS % % % % % 
    elseif strcmpi(algo,'SVDLS')
        rsvd = 'full';
        [i,j] = CUR_LeverageScore(target.A, rank, rsvd, target.U, target.V);
    elseif strcmpi(algo, 'RSVDLS')
        rsvd = 'gauss';
        [i,j] = CUR_LeverageScore(target.A, rank, rsvd);
    elseif strcmpi(algo, 'CSSVDLS')
        randmat = 'sparse3';
        [i,j] = CUR_LeverageScore(target.A, rank, randmat);
    elseif strcmpi(algo, 'RSVDLSstream')
        randmat = 'gauss';
        stream = 1;
        [i,j] = CUR_LeverageScore(target.A, rank, randmat, stream);
% % % % % default % % % % %     
    else % CPQR
        randmat = 'gauss';
        [i,j] = CUR_ID(target.A, rank, randmat);
    end
end
              