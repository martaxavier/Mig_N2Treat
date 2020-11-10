% ------------------------------------------------------------
% Go through metrics
% ------------------------------------------------------------

for m = 1 : length(metrics)
    
    metric = metrics(metric);
    get_metric_pars;
   
    % Continue if metric is not yet supported for model fitting 
    if max(contains(metrics_not_yet_supported,eeg_metric))
        continue
    end
   
    % Load estimated order of the ACF
    if strcmp(bold_shift,'deconv')
        acf_order_in = 'acf_order_deconv.mat';
    else
        acf_order_in = 'acf_order.mat';
    end
    load(fullfile(path_pars_in, ...
        acf_order_in),'acf_order');
    
    % Compute average ACF order
    order_avg = mean(acf_order.order);
    
    % ------------------------------------------------------------
    % Go through subjects and get X and Y
    % ------------------------------------------------------------
    for s = 1 : length(subjects)

        subject = subjects(s);

        display(strcat('Estimating CV pars for'," ", ...
            subject,','," ",metric,','," ",cv,'...'));

        % Define input EEG and BOLD data, according
        % to current metric
        eeg_in = strcat(eeg_metric,'_', ...
            'eeg_feature','_',eeg_shift,'.txt');
        bold_in = strcat('bold_processed', ...
            '_',bold_shift,'.txt');

        if contains(eeg_in,'_.')
            eeg_in = replace(eeg_in,'_.','.');
        end
        if contains(bold_in,'_.')
            bold_in = replace(bold_in,'_.','.');
        end

        % Read subject eeg and bold files 
        X(:,:,s) = dlmread(strcat(path_eeg_in,eeg_in));
        Y(:,:,s) = dlmread(strcat(path_bold_in,bold_in));

    end

    siz_X = [size(X,1) prod(dim(1:2)) dim(3)];

    % ------------------------------------------------------------
    % Get stats and best k-v pair for all subjects 
    % ------------------------------------------------------------

    [bestpair,stats] = get_cv_pars(X,Y,'cv',cv_method,...
        'regress',reg_model,'autocorr',order_avg,'sizx',siz_X);

    % ------------------------------------------------------------
    % Save output matrices with optimal parameters  
    % ------------------------------------------------------------

    % Save best K-V pair in outputstruct
    cvpars.K = bestpair(1); cvpars.V = bestpair(2);

    % Save matrix of best parameters for this cv method  
    cv_pars_out = strcat(cv_method,'_',metric,'.mat');
    save(fullfile(path_pars_out,cv_pars_out),'cvpars');    

end

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

function [bestpair,stats] = get_cv_pars(X,Y,varargin)

%   Input parameters:
%
%     X                  numeric matrix, NxPxS
%     Y                  numeric vector of length NxS
%
%   Optional input parameters:  
%
%     'cv'             The cv method to be tested 
%                      ('regular','nondep','blocked')
%     'regress'        The regression method ('elasticnet','l21_1')
%     'autocorr'       The order of the ACF
%     'sizx'           The size of the problem
%
%   Outputs:
%
%       bestpairs       optimal K-V pair. K is the total number of
%                       folds, V is the validation to learn set ratio
%
%
%       stats            struct that contains information about the
%                        sequence of model fits corresponding to the
%                        the searched parameters. stats contains the
%                        following fields:
%  
%       'nmse'           normalized mean squared error of the 
%                        fits for each of the combinations tested 
%       'df'             degrees of freedom of the fits for 
%                        each of the combinations tested 
%       'bic'            bayesian information criterion of the 
%                        fits for each of the combinations tested 

tic 

% ------------------------------------------------------------
% Sanity check and process input parameters 
% ------------------------------------------------------------

% X a real 3D matrix
if  length(size(X)) ~= 3 || ~isreal(X)
    error('X is not a real 3D matrix');
end

% Y is a real 3D matrix
if  length(size(Y)) ~= 3 || ~isreal(Y)
    error('Y is not a real 3D matrix');
end

% Make sure dim 2 is a singleton
if size(Y,2) ~= 1
    Y = permute(Y,[2 1 3]);
end

% Sizes of X and Y agree
if size(Y,1) ~= size(X,1) || size(Y,3) ~= size(X,3)
    error('X and Y dimensions do not agree');
end

% Define # time-points, # features, # subjects
[n_pnts,~,n_subjs] = size(X);

% ------------------------------------------------------------
% Sanity check and process optional parameters 
% ------------------------------------------------------------

% varargin is a cell that contains the input arguments
% in string format 
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% Assign default values for each optional parameter
pnames = {'cv' 'regress' 'autocorr' 'sizx'};
dflts  = {'regular' 'elasticnet' [n_pnts 31*6 4]};

% Assign variables corresponding to optional parameters 
[cv, reg_model, order, siz_X] = internal.stats.parseArgs...
    (pnames, dflts, varargin{:});

% ------------------------------------------------------------
% Choose the K and V values to be tested 
% ------------------------------------------------------------

% Typical values of K
% range from 5 to 15 
K = 10:5:15; 
n_k = length(K);

% Typical values of V 
% range from 0.1 to 0.3
V = 0.1:0.1:0.3; 
n_v = length(V);

% Linear indices 
idxs = 1:n_k*n_v;
n_pairs = length(idxs);

% Subscript indices 
[k_idxs, v_idxs] = ...
    ind2sub([n_k,n_v],idxs);

nmse = zeros(n_subjs,n_pairs); 
df = nmse; bic = nmse;

% ------------------------------------------------------------
% Go through K-V pairs 
% ------------------------------------------------------------
for ii = 1 : n_pairs
    
    k = K(k_idxs(ii));
    v = V(v_idxs(ii));
        
    % ------------------------------------------------------------
    % Go through  subjects 
    % ------------------------------------------------------------
    for s = 1 : n_subjs
        
        EEG = X(:,:,s);
        BOLD = Y(:,:,s);
        
        % Center and scale EEG
        [EEG,~,~] = zscore(EEG,1);

        % Center BOLD
        muBOLD = mean(BOLD);
        BOLD = bsxfun(@minus,BOLD,muBOLD);
        
        switch cv 
            case 'regular'
                [model,~] = kfold_cv_par_v2(EEG,BOLD, ...
                    'k',k,'val2learn',v,'regress',reg_model);  
            case 'nondep'
                [model,~] = kfold_cv_nondep_par_v3(EEG,BOLD,'k',k,...
                    'val2learn',v,'regress',reg_model,'autocorr', ...
                    order,'sizx',siz_X);     
                    
            case 'blocked'
                [model,~] = kfold_cv_blocked_par_v2(EEG,BOLD,'k',k,...
                    'val2learn',v,'regress',reg_model,'autocorr', ...
                    order,'sizx',siz_X);     
        end 
    
        nmse(s,ii) = model.nmse;
        df(s,ii) = model.df;
        bic(s,ii) = model.bic;
        
        
    end % finish looping through subjects 

end % finish looping through k-v pairs 

% ------------------------------------------------------------
% Interpolate remaining K-V pairs
% ------------------------------------------------------------

% Average model results across subject dimension
nmse_avg = mean(nmse,1);
df_avg = mean(df,1); 
bic_avg = mean(bic,1);

% Reshape results into matrices for plotting 
nmse_mat = reshape(nmse_avg,[n_k,n_v]);
df_mat = reshape(df_avg,[n_k,n_v]);
bic_mat = reshape(bic_avg,[n_k,n_v]);

% Plot surf of nmse, df and bic 
[K_2d,V_2d] = meshgrid(K,V);

figure;surf(K_2d,V_2d,nmse_mat);
title('NMSE'); xlabel('K'); ylabel('V');

figure;surf(K_2d,V_2d,df_mat);
title('DF'); xlabel('K'); ylabel('V');

figure;surf(K_2d,V_2d,bic_mat);
title('BIC'); xlabel('K'); ylabel('V');

% ------------------------------------------------------------
% Find best K-V pair 
% ------------------------------------------------------------

% Sort indices according to bic 
[~,sortedindices] = sort(bic_avg);
bestpair = [K(k_idxs((sortedindices(1)))),...
    V(v_idxs(sortedindices(1)))];


% ------------------------------------------------------------
% Prepare outputs 
% ------------------------------------------------------------

% Prepare stats struct 
stats.k = K(k_idxs((sortedindices)))';
stats.v = V(v_idxs((sortedindices)))';
stats.nmse = nmse_avg(sortedindices)';
stats.df = df_avg(sortedindices)';
stats.bic = bic_avg(sortedindices)';
stats.summary = struct2table(stats);

toc

end 