function [LL_all,CloneHMM_paras,p_states,aCN,segments] = CloneHMM_screening(init_CloneHMM_paras,thres1,max_iter1,verbose)

global data_lrc_sep

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_CloneHMM_paras{1};
transmat_all = init_CloneHMM_paras{2};
o_all = init_CloneHMM_paras{3};
sigma_all = init_CloneHMM_paras{4};
indivec_all = init_CloneHMM_paras{5};

LL_all = [];
CloneHMM_paras = cell(1,5); 
if nargout > 2
    p_states = [];
    aCN = zeros(size(data_lrc_sep{1},1),length(o_all));
    segments = cell(size(data_lrc_sep{1},1),length(o_all));
end

for i = 1:length(o_all)
    %1x1 cell
    init_CloneHMM_paras(1) = prior_all(i);
    init_CloneHMM_paras(2) = transmat_all(i);
    init_CloneHMM_paras(3) = o_all(i);
    init_CloneHMM_paras(4) = sigma_all(i);
    init_CloneHMM_paras(5) = indivec_all(i);
    
    [LL,prior,transmat,o,sigma,iterations] = CloneHMM_estimate_paras(init_CloneHMM_paras,thres1,max_iter1,verbose);
        
    LL_all = [LL_all LL(end)];
    CloneHMM_paras{1} = [CloneHMM_paras{1} {prior}];
    CloneHMM_paras{2} = [CloneHMM_paras{2} {transmat}];
    CloneHMM_paras{3} = [CloneHMM_paras{3} {o}];
    CloneHMM_paras{4} = [CloneHMM_paras{4} {sigma}];
    CloneHMM_paras{5} = [CloneHMM_paras{5} init_CloneHMM_paras(5)];
    
    if nargout > 2
        [temp,aCN(:,i),segments(:,i)] = CloneHMM_process_results();
        p_states = [p_states {temp}];
    end

    if verbose
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, iterations:' num2str(iterations)]);
        disp(['sigma:' num2str(sigma) ', o:' num2str(o) ', LL:' num2str(LL(end),'%5.1f')]);
        disp('--------------- screening report -----------------')
    end
    
end