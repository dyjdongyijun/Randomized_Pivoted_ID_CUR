function test_CUR_problem_size(rank, algos, tag, target, varargin)
% algos = cell of strings representing CUR algorithms being tested on
% target = 'gpsd': {prob_sizes, sig}
% target = otherwise: {nvertex, nedge}
%% Construct target matrix
    if nargin < 1
        rank = 20;
    end
    if nargin < 2 || isempty(algos)
        algos = {'DetCPQR',...
                'CPQR',...
                'CSCPQR',...
                'CPQRstream',...
                ...
                'DetLUPP',...
                'LUPP',...
                'CSLUPP',...
                'LUPPstream',...
                ...
                'SVDDEIM',...
                'RSVDDEIM',...
                'CSSVDDEIM',...
                'RSVDDEIMstream',...
                ...
                'LUCP',...
                'ACA',...(stream)
                'CSLUCP',...(stream)
                ...
                'SVDLS',...
                'RSVDLS',...
                'CSSVDLS',...
                'RSVDLSstream'};
    end
    
    if nargin < 3, tag = []; end
    
    if nargin < 4 % default: dense GPSD
        prob_sizes = 400:400:8000;
        m_max = 8000;
        amp = 2; 
        len = 10;
        sig = 1./(1:m_max); 
        sig(1:len) = sig(1:len)*amp;
        
        test_CUR_prob_size_GPSD(algos, prob_sizes, rank, sig);
    
    elseif strcmpi(target, 'gpsd')
        if ~isempty(varargin)
            prob_sizes = varargin{1};
        else
            prob_sizes = 400:400:8000;
        end
        if length(varargin) >= 2
            sig = varargin{2};
        else
            m_max = prob_sizes(end);
            amp = 2; 
            len = 10;
            sig = 1./(1:m_max); 
            sig(1:len) = sig(1:len)*amp;
        end
           
        test_CUR_prob_size_GPSD(algos, prob_sizes, rank, sig, tag);
    
    else %if strcmpi(target,'laplacian')
        if ~isempty(varargin)
            nvertex = varargin{1};
        else
            nvertex = 1000:1000:10000;
        end
        if length(varargin) >= 2
            nedge = varargin{2};
        else
            nedge = 5*nvertex;
        end
        
        test_CUR_prob_size_Laplacian(algos, nvertex, nedge, rank, tag);
    end
end

%% Test on dense GPSD matrix
function test_CUR_prob_size_GPSD(algos, prob_sizes, rank, sig, tag)
    k = rank;
    parameters = struct('description','dense Gaussian SPD',....
                    'rank',k,...
                    'problem_sizes',prob_sizes,...
                    'sigma',sig);
    time = struct();
    for idx = 1:length(algos)
        algo = algos{idx};
        time.(algo) = size(prob_sizes);
    end
    %% Latency
    m = 50;
    target = TargetMatGenerator('gspd',m,m,sig);
    for idx = 1:length(algos)
        algo = algos{idx};
        tic;
        [~,~] = str2curalgo(algo, target, k);
        toc;
    end
    
    %% Experiment
    for t = 1:length(prob_sizes)
        m = prob_sizes(t);
        target = TargetMatGenerator('gspd',m,m,sig);

        for idx = 1:length(algos)
            algo = algos{idx};
            tic;
            [~,~] = str2curalgo(algo, target, k);
            time.(algo)(t) = toc;
        end
        
        fprintf('m = %d / %d finished \n', m, prob_sizes(end))
    end

    %% save data
    if isempty(tag)
        tag = sprintf('k%d_dense-gspd_n%d', rank, prob_sizes(end));
    end
    probsizes_file = sprintf('probsizes_%s.mat', tag);
    time_file = sprintf('time_%s.mat', tag);
    save(probsizes_file,'-struct','parameters')
    save(time_file,'-struct','time')
    fprintf('write out: %s %s \n', probsizes_file, time_file)
end

%% Test on sparse Laplacian
function test_CUR_prob_size_Laplacian(algos, nvertex, nedge, rank, tag)
    k = rank;
    ns = nvertex;
    ms = nedge;
    parameters = struct('description','Random uniformly weighted Laplacian',....
                        'rank',k,...
                        'num_V',ns,...
                        'num_E',ms);
    time = struct();
    for idx = 1:length(algos)
        algo = algos{idx};
        time.(algo) = zeros(size(ns));
    end
    %% Latency
    n = 2*k; m = 5*n;
    target = TargetMatGenerator('laplacian',n,m);
    for idx = 1:length(algos)
        algo = algos{idx};
        tic;
        [~,~] = str2curalgo(algo, target, k);
        toc;
    end
    
    %% Experiment
    for t = 1:length(ns)
        n = ns(t); m = ms(t);
        target = TargetMatGenerator_raw('laplacian',n,m);
        
        for idx = 1:length(algos)
            algo = algos{idx};
            tic;
            [~,~] = str2curalgo(algo, target, k);
            time.(algo)(t) = toc;
        end
        
        fprintf('%d / %d finished \n', t, length(ns))
    end

    %% save data
    if isempty(tag)
        tag = sprintf('k%d_sparse-laplacian_n%d', rank, ns(end));
    end
    probsizes_file = sprintf('probsizes_%s.mat', tag);
    time_file = sprintf('time_%s.mat', tag);
    save(probsizes_file,'-struct','parameters')
    save(time_file,'-struct','time')
    fprintf('write out: %s %s \n', probsizes_file, time_file)
end

%% % % % % % % % % % % % test functions for each algo % % % % % % % % % 
function [i,j,U] = str2curalgo(algo, target, rank)
% % % % % LUPP % % % % % 
    if strcmpi(algo, 'DetLUPP')
        sketch = [];
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'LUPP')
        sketch = 'gauss';
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'CSLUPP')
        sketch = 'sparse3';
        [i,j] = CUR_LUPP(target.A, rank, sketch);
    elseif strcmpi(algo, 'LUPPstream')
        sketch = 'gauss';
        stream = 1;
        [i,j] = CUR_LUPP(target.A, rank, sketch, stream);
% % % % % CPQR % % % % % 
    elseif strcmpi(algo, 'DetCPQR')
        rsvd = 0;
        [i,j] = CUR_ID(target.A, rank, rsvd);
    elseif strcmpi(algo, 'CPQR')
        randmat = 'gauss';
        [i,j] = CUR_ID(target.A, rank, randmat);
    elseif strcmpi(algo, 'CSCPQR')
        randmat = 'sparse3';
        [i,j] = CUR_ID(target.A, rank, randmat);
    elseif strcmpi(algo, 'CPQRstream')
        randmat = 'gauss';
        [i,j] = CUR_ID_streaming(target.A, rank, randmat);
    elseif strcmpi(algo, 'CPQRstreamCUR')
        randmat = 'gauss';
        [i,j,U] = CUR_ID_streaming(target.A, rank, randmat);
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
