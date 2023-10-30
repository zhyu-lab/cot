function SCHMM_paras = SCHMM_main(init_SCHMM_paras,thres_EM,max_iter,verbose)

%--------------------------- screening -------------------------
global candi_cns
global clamp_thres
global mc_w
global NoSolutionFlag
clamp_thres = 1-1e-3;
mc_w = 0.8;

%initialize parameters
thres_del = 0.009;

init_SCHMM_paras = SCHMM_init_paras(init_SCHMM_paras,0);
[LL,SCHMM_paras,p_states,aCN] = SCHMM_screening...
    (init_SCHMM_paras,thres_EM,max_iter,verbose);

tv_del = candi_cns < 1;
tv_valid = zeros(1,length(p_states))==1;
for i = 1:length(p_states)
    tv = sum(p_states{i}(:,tv_del),2) < thres_del & aCN(:,i) < 4.5;
    if sum(tv) == length(tv)
        tv_valid(i) = 1;
    end
end

NoSolutionFlag = false;
if ~any(tv_valid)
    tv_valid = ~tv_valid;
    warning('Can not find a feasible solution with pre-defined criteria!');
    NoSolutionFlag = true;
end
candi_indx = find(tv_valid);
candi_ll = LL(tv_valid);
candi_acn = mean(aCN(:,tv_valid),1);

[temp,I] = max(candi_ll);
best_indx = candi_indx(I);

[temp,indxs] = sort(candi_ll,'descend');
scores = 1;
acns = candi_acn(indxs(1));
for i = 2:length(candi_acn)
    j = indxs(i-1);
    pre_indxs = indxs(1:i-1);
    k = indxs(i);
    if abs(candi_ll(j)-candi_ll(k)) <= 2 && abs(candi_acn(j)-candi_acn(k)) <= 0.1
        continue;
    end
    score = (candi_acn(pre_indxs)-candi_acn(k)+eps)./(candi_ll(pre_indxs)-candi_ll(k)+eps);
    scores = [scores score(end)];
    acns = [acns candi_acn(k)];
    if sum(score >= 0.02) == i-1
        best_indx = candi_indx(k);
    end
end

% tmp = (aCN(best_indx)-aCN(tv)+eps)./(LL(best_indx)-LL(tv)+eps);
%disp(['acn: ' num2str(acns)])
%disp(['score: ' num2str(scores)])
%disp(['ll: ' num2str(LL)])

% Now, use the optimal parameters to call CNA
%-------------------------------------------------------------------
init_SCHMM_paras = SCHMM_init_paras(SCHMM_paras,best_indx);
[temp,SCHMM_paras] = SCHMM_screening(init_SCHMM_paras,5*thres_EM,20,verbose);
%-------------------------------------------------------------------


function SCHMM_paras = SCHMM_init_paras(init_SCHMM_paras,best_indx)
%this function is used to initialize/process parameters for HMM training,
%best_indx is used to indicate which parameter configuration (ususally
%multiple generated in previous screening procedure) are selected. If
%best_indx equals 0, parameters will be initialized

global clamp_thres
global var_l
global candi_cns

S = length(candi_cns);
SCHMM_paras = cell(1,5);
%parameter initialization
if best_indx == 0
    %---w---
    if isempty(init_SCHMM_paras{3})
        o_0 = [-1.2 -1 -0.6 -0.3 0];
    else
        o_0 = init_SCHMM_paras{3};
    end
    
    N = length(o_0);
    SCHMM_paras{3} = mat2cell(o_0,1,ones(1,N));
    
    %---sigma--- 
    if isempty(init_SCHMM_paras{4})
        sigma_0 = sqrt(var_l);
    else
        sigma_0 = init_SCHMM_paras{4};
    end
    SCHMM_paras{4} = repmat({sigma_0},1,N);

    %---pi---
    if isempty(init_SCHMM_paras{1})
        prior_0 = 1/(S)*ones(S,1);
    else
        prior_0 = init_SCHMM_paras{1};
    end
    SCHMM_paras{1} = repmat({prior_0},1,N);

    %---A---
    if isempty(init_SCHMM_paras{2})
        transmat_0 = norm_trans(ones(S,S),clamp_thres);
    else
        transmat_0 = init_SCHMM_paras{2};
    end
    SCHMM_paras{2} = repmat({transmat_0},1,N);
    
    %---indicator vector---
    if isempty(init_SCHMM_paras{5}) %indicator vector: '1' for update '0' fixed
        adj_all = ones(1,4);
    else
        adj_all = init_SCHMM_paras{5};
    end
%     adj_all(4) = 0;
    SCHMM_paras{5} = repmat({adj_all},1,N);
else %parse the results from previous screening
    for i = 1:length(best_indx)
        %--pi--
        SCHMM_paras{1} = [SCHMM_paras{1} init_SCHMM_paras{1}(best_indx(i))];
        %--A--
        SCHMM_paras{2} = [SCHMM_paras{2} init_SCHMM_paras{2}(best_indx(i))];
         %--o--
        SCHMM_paras{3} = [SCHMM_paras{3} init_SCHMM_paras{3}(best_indx(i))];
        %--sigma--
        SCHMM_paras{4} = [SCHMM_paras{4} init_SCHMM_paras{4}(best_indx(i))];
        %--indicator vector--
        SCHMM_paras{5} = [SCHMM_paras{5} init_SCHMM_paras{5}(best_indx(i))];
    end
end
