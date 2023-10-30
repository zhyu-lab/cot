function [LL_all,SCHMM_paras,p_states,aCN,segments] = SCHMM_screening(init_SCHMM_paras,thres1,max_iter1,verbose)

global data_lrc_sep

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_SCHMM_paras{1};
transmat_all = init_SCHMM_paras{2};
o_all = init_SCHMM_paras{3};
sigma_all = init_SCHMM_paras{4};
indivec_all = init_SCHMM_paras{5};

LL_all = [];
SCHMM_paras = cell(1,5); 
if nargout > 2
    p_states = [];
    aCN = zeros(size(data_lrc_sep{1},1),length(o_all));
    segments = cell(size(data_lrc_sep{1},1),length(o_all));
end

for i = 1:length(o_all)
    %1x1 cell
    init_SCHMM_paras(1) = prior_all(i);
    init_SCHMM_paras(2) = transmat_all(i);
    init_SCHMM_paras(3) = o_all(i);
    init_SCHMM_paras(4) = sigma_all(i);
    init_SCHMM_paras(5) = indivec_all(i);
    
    [LL,prior,transmat,o,sigma,iterations] = SCHMM_estimate_paras(init_SCHMM_paras,thres1,max_iter1,verbose);
        
    LL_all = [LL_all LL(end)];
    SCHMM_paras{1} = [SCHMM_paras{1} {prior}];
    SCHMM_paras{2} = [SCHMM_paras{2} {transmat}];
    SCHMM_paras{3} = [SCHMM_paras{3} {o}];
    SCHMM_paras{4} = [SCHMM_paras{4} {sigma}];
    SCHMM_paras{5} = [SCHMM_paras{5} init_SCHMM_paras(5)];
    
    if nargout > 2
        [temp,aCN(:,i),segments(:,i)] = SCHMM_process_results();
        p_states = [p_states {temp}];
    end

    if verbose
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, iterations:' num2str(iterations)]);
        disp(['sigma:' num2str(sigma) ', o:' num2str(o) ', LL:' num2str(LL(end),'%5.1f')]);
        disp('--------------- screening report -----------------')
    end
    
end