function [p_states,aCN,segments_all] = CloneHMM_process_results()

%-----------------------------------------------------
%------overall information of the clone------
%p_states: proportions of all hidden states
%aCN: averaged copy number
%segments: copy number segmentation results

global candi_cns
global gamma_sep

%initialize intermediate variables
exp_num_states = [];
pos_dist = [];

num_cell = length(gamma_sep{1});
segments_all = cell(num_cell,1);
pos_dist_all = cell(num_cell,1);
for c = 1:num_cell
    segments = [];
    pos_dist = [];
    for ex = 1:length(gamma_sep) %for the ith chromosome
        post_probs = gamma_sep{ex}{c};

        %---handle p_states and num_loci---
        if isempty(exp_num_states) %initialization
            exp_num_states = zeros(size(post_probs,1),1);
        end
        exp_num_states = exp_num_states+sum(post_probs,2);

        %---handle MAP states---
        %output predicted MAP states
        [temp,MAP_state] = max(post_probs,[],1);

        results = CloneHMM_segment_results(MAP_state);
        segments = [segments; ones(size(results,1),1)*ex results];    
        pos_dist = [pos_dist; (results(:,2)-results(:,1)+1)]; 
    end
    segments_all{c} = segments;
    pos_dist_all{c} = pos_dist;
end

%---handle p_states---
p_states = zeros(num_cell,length(candi_cns));
for c = 1:num_cell
    segments = segments_all{c};
    pos_dist = pos_dist_all{c};
    for i = 1:length(candi_cns)
        tv = segments(:,4) == i;
        if sum(tv) > 0
            p_states(c,i) = sum(pos_dist(tv))/sum(pos_dist);
        end
    end
end

%---handle aCN---
aCN = p_states*candi_cns';

end