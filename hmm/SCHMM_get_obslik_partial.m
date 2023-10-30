function [obslik,condi_probs_fluct] = SCHMM_get_obslik_partial(data_lrc,o,sigma)

global candi_cns

N = length(data_lrc); %number of data points

mu_l = log2(candi_cns/2)+o;

S = length(mu_l);
obslik = zeros(S,N);
condi_probs_fluct = zeros(S,N);

fluct_prob = 1e-5;

for i = 1:length(mu_l)
    obslik_lrc = SCHMM_eval_pdf_lrc(data_lrc,mu_l(i),sigma);
    if candi_cns(i) == 2
        obslik_lrc = 1.06*obslik_lrc;
    end
    if candi_cns(i) == 1
        obslik_lrc = 1.2*obslik_lrc;
    end
    obslik(i,:) = (1-fluct_prob)*obslik_lrc+fluct_prob/6;
    condi_probs_fluct(i,:) = (fluct_prob/6)./obslik(i,:);
end

end
