function SCHMM(lrcFile,labelFile,outputDir,maxCN)
% 02/05/2023

%------Input and output------%
% lrcFile: file containing log read counts of all cells
% labelFile: file containing cluster labels of all cells
% outputDir: directory of output files
% maxCN: maximum copy number to consider

global current_version
current_version = '1.0';

if nargin < 3
    error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m'] );
end
if nargin < 4
	maxCN = 10;
end

if strcmp(class(maxCN),'char') == 1
    maxCN = str2double(maxCN);
end

disp(['maximum copy number is set to ' num2str(maxCN)]);

thres_EM = 1e-5;
max_iter = 50;
verbose = 1;
                                                                   
%initialization of global variable
global data_lrc_all;
global data_lrc_sep
global data_bin_sep
global bin_size
global var_l

global candi_cns

candi_cns = [0.001 1:maxCN];

% record time
tic

mkdir(outputDir);

disp(['SCHMM (version ' current_version ') is loading...'])

%--------------load data--------------------
[data_lrc_all,data_chr_all,data_bin_all,bin_size,cell_labels] = SCHMM_load_data(lrcFile,labelFile);
[num_cell,num_bin] = size(data_lrc_all);
cluster_ids = unique(cell_labels);
num_cluster = length(cluster_ids);

chromosomes = reshape(unique(data_chr_all),1,[]);
num_chrs = length(chromosomes);

disp('----------Call CNAs now----------');

data_paras_all = zeros(num_cluster,2);
data_cn_all = zeros(num_cell,num_bin);
data_acn_all = zeros(num_cell,1);
data_segments_all = cell(num_cell,1);
for c = 1:num_cluster
%     if cluster_ids(c) ~= 5
%         continue;
%     end
    tv = cell_labels == cluster_ids(c);
    cell_indxs = find(tv);
    data_lrc_cluster = data_lrc_all(tv,:);
    var_l = var(data_lrc_cluster(:));

    %-------divide into different chromosomes----------
    data_bin_sep = cell(num_chrs,1);
    data_lrc_sep = cell(num_chrs,1);

    for i = 1:num_chrs
        tv = ismember(data_chr_all,chromosomes(i));
        data_bin_sep{i} = data_bin_all(tv);
        data_lrc_sep{i} = data_lrc_cluster(:,tv);
    end
    clear data_lrc_cluster;
    
    init_SCHMM_paras = [{[]},{[]},{[]},{[]},{[]}]; % initial parameters:pie,transmat,sigma,o
    %------------------ call HMM --------------------
    disp(['----------Call CNAs for cluster ' num2str(c) ' ----------']);
    SCHMM_paras = SCHMM_main(init_SCHMM_paras,thres_EM,max_iter,verbose);
    o = SCHMM_paras{3}{1};
    sigma = SCHMM_paras{4}{1};
    data_paras_all(c,:) = [o sigma];
 
    [p_states,aCN,segments_all] = SCHMM_process_results();

    for s = 1:length(segments_all)
        segments = segments_all{s};
        cn_segs_all = zeros(size(segments,1),1); % copy number
        bp_len = 0;
        cn_w = 0;
        j = 1;

        data_segs_cell = zeros(size(segments,1),4);
        for i = 1:size(segments,1)
            chr_indx = segments(i,1);
            s_indx = segments(i,2);
            e_indx = segments(i,3);
            state_indx = segments(i,4);
            St_pos = (data_bin_sep{chr_indx}(s_indx)-1)*bin_size+1;
            Ed_pos = data_bin_sep{chr_indx}(e_indx)*bin_size;

            data_lrc = data_lrc_sep{chr_indx}(s,s_indx:e_indx);

            lrc_mean = median(data_lrc);

            CN = round(2^(lrc_mean-o+1));

            % assign copy number to each segment
            if CN ~= candi_cns(state_indx)
                cn_segs_all(i) = CN;
            else
                cn_segs_all(i) = candi_cns(state_indx);
            end

            bp_len = bp_len+(Ed_pos-St_pos+1);   
            cn_w = cn_w+(Ed_pos-St_pos+1)*cn_segs_all(i);

            data_cn_all(cell_indxs(s),j:j+e_indx-s_indx) = cn_segs_all(i);
            j = j+e_indx-s_indx+1;

            data_segs_cell(i,:) = [chromosomes(chr_indx) St_pos Ed_pos cn_segs_all(i)];

        end

        data_acn_all(cell_indxs(s)) = cn_w/bp_len;
        data_segments_all{cell_indxs(s)} = data_segs_cell;
    end
end

%---save results---%
o_fid = fopen([outputDir '/paras.csv'],'w');
if o_fid == -1
    error(['Can not open file ' outputDir '/paras.csv for saving results.']);
end
fprintf(o_fid,'baseline shift,sigma\n');
for k = 1:num_cluster
    fprintf(o_fid,'%.6f,%.6f\n',data_paras_all(k,1),data_paras_all(k,2));
end
fclose(o_fid);

o_fid = fopen([outputDir '/segments.csv'],'w');
if o_fid == -1
    error(['Can not open file ' outputDir '/segments.csv for saving results.']);
end
fprintf(o_fid,'Cell,Chr,StartPos,EndPos,CN\n');
for k = 1:num_cell
    segments = data_segments_all{k};
    for i = 1:size(segments,1)
        chr = segments(i,1);
        s_pos = segments(i,2);
        e_pos = segments(i,3);
        cn = segments(i,4);
        fprintf(o_fid,'%d,%d,%d,%d,%d\n',k,chr,s_pos,e_pos,cn);
    end
end
fclose(o_fid);

o_fid = fopen([outputDir '/copynumber.csv'],'w');
if o_fid == -1
    error(['Can not open file ' outputDir '/copynumber.csv for saving results.']);
end
for i = 1:num_bin-1
    fprintf(o_fid,'%d:%d-%d,',data_chr_all(i),(data_bin_all(i)-1)*bin_size+1,data_bin_all(i)*bin_size);
end
fprintf(o_fid,'%d:%d-%d\n',data_chr_all(end),(data_bin_all(end)-1)*bin_size+1,data_bin_all(end)*bin_size);
for k = 1:num_cell
    for i = 1:num_bin-1
        fprintf(o_fid,'%d,',data_cn_all(k,i));
    end
    fprintf(o_fid,'%d\n',data_cn_all(k,end));
end
fclose(o_fid);

o_fid = fopen([outputDir '/ploidy.csv'],'w');
if o_fid == -1
    error(['Can not open file ' outputDir '/ploidy.csv for saving results.']);
end
for k = 1:num_cell-1
    fprintf(o_fid,'%.3f,',data_acn_all(k));
end
fprintf(o_fid,'%.3f\n',data_acn_all(end));
fclose(o_fid);

t_all = toc;
disp(['-----Finish calling CNAs, totally ' num2str(t_all/60) ' minites were used-----'])

end