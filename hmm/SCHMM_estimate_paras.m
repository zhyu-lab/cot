function [LL,prior,transmat,o,sigma,nrIterations] = ...
    SCHMM_estimate_paras(init_SCHMM_paras,thresh,max_iter,verbose)

global clamp_thres

previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];

prior = init_SCHMM_paras{1};
transmat = init_SCHMM_paras{2};
o = init_SCHMM_paras{3};
sigma = init_SCHMM_paras{4};

while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [loglik,exp_num_trans,exp_num_visits1,o_u,sigma_u] = ...
        SCHMM_compute_ess(prior,transmat,o,sigma);
    
    converged = em_converged_m(loglik,previous_loglik,verbose,thresh);
    
    % update parameters
    if init_SCHMM_paras{5}(1)
        prior = norm_trans(exp_num_visits1',0)';
    end
    if init_SCHMM_paras{5}(2) && ~isempty(exp_num_trans)
        % clamp_thres = 1-1e-4;
        transmat = norm_trans(exp_num_trans,clamp_thres);
    end
    if init_SCHMM_paras{5}(3) %update o here
        o = o_u;
    end
    if init_SCHMM_paras{5}(4) %update sigma here
        sigma = sigma_u;
    end
    
    if verbose
        disp(['sigma:' num2str(sigma) ', o:' num2str(o)]);
        fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik);
    end
    
    num_iter =  num_iter + 1;
    previous_loglik = loglik;
    LL = [LL loglik];
end
nrIterations = num_iter - 1;

end

%--------------------------------------------------------------------------
function [loglik,exp_num_trans,exp_num_visits1,o_u,sigma_u] = ...
    SCHMM_compute_ess(prior,transmat,o,sigma)

global data_lrc_sep
global gamma_sep
global condi_probs_fluct_sep

num_chrs = length(data_lrc_sep);
S_all = size(transmat,1); % number of states 
exp_num_trans = zeros(S_all,S_all);
exp_num_visits1 = zeros(S_all,1);

%-----------------------E step-----------------------------
gamma_sep = cell(1,num_chrs);
condi_probs_fluct_sep = cell(1,num_chrs);
loglik = 0;

for i = 1:num_chrs
    data_lrc_chr = data_lrc_sep{i};
    num_cell = size(data_lrc_chr,1);
    for c = 1:num_cell
        % conditional probabilities
        [obslik,condi_probs_fluct] = SCHMM_get_obslik(data_lrc_chr(c,:),o,sigma);
        % Forward and Backward algorithm
        [alpha,gamma,current_ll,beta,xi_summed] = Forward_Backward_Algorithm(prior,transmat,obslik);
        clear alpha beta;
        
        loglik = loglik + current_ll;
        exp_num_trans = exp_num_trans + xi_summed;
        exp_num_visits1 = exp_num_visits1 + gamma(:,1);

        gamma_sep{i} = [gamma_sep{i} {gamma}];
        condi_probs_fluct_sep{i} = [condi_probs_fluct_sep{i} {condi_probs_fluct}];
        clear gamma condi_probs_fluct;
    end
end

%-----------------------M step-----------------------------
%update sigma
sigma_u = SCHMM_update_sigma(o,sigma);

%update o
o_u = SCHMM_update_o(o);

end

%--------------------------------------------------------------------------
function sigma_u = SCHMM_update_sigma(o,sigma)

global candi_cns
global data_lrc_sep
global gamma_sep
global condi_probs_fluct_sep

num_chrs = length(data_lrc_sep); % each row is a sample
mu_l = log2(candi_cns/2)+o;

numerator = 0;
denominator = 0;

for ex = 1:num_chrs
    data_lrc_chr = data_lrc_sep{ex};
    num_cell = size(data_lrc_chr,1);
    for c = 1:num_cell
        gamma = gamma_sep{ex}{c}.*(1-condi_probs_fluct_sep{ex}{c});
        for i = 1:length(mu_l)
            numerator = numerator+gamma(i,:)*((data_lrc_chr(c,:)-mu_l(i)).^2)';
            denominator = denominator+sum(gamma(i,:)); 
        end
    end
end

sigma_u = sqrt(numerator/denominator);
if isnan(sigma_u)
    sigma_u = sigma;
end

end

%--------------------------------------------------------------------------
function o_u = SCHMM_update_o(o)

global candi_cns
global data_lrc_sep
global gamma_sep
global condi_probs_fluct_sep

num_chrs = length(data_lrc_sep); % each row is a sample

temp = log2(candi_cns/2);

numerator = 0;
denominator = 0;

for ex = 1:num_chrs
    data_lrc_chr = data_lrc_sep{ex};
    num_cell = size(data_lrc_chr,1);
    for c = 1:num_cell
        gamma = gamma_sep{ex}{c}.*(1-condi_probs_fluct_sep{ex}{c});
        for i = 1:length(candi_cns)
            numerator = numerator+gamma(i,:)*(data_lrc_chr(c,:)-temp(i))';
            denominator = denominator+sum(gamma(i,:)); 
        end
    end
end

o_u = numerator/denominator;
if isnan(o_u)
    o_u = o;
end

end

