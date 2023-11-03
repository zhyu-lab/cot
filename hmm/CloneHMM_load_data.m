function [data_lrc_all,data_chr_all,data_bin_all,bin_size,cell_labels] = CloneHMM_load_data(lrcFile,labelFile)

fid = fopen(labelFile,'r');
fgetl(fid);
line = fgetl(fid);
fields = regexp(line,',','split');
cell_labels = str2double(fields);
fclose(fid);

fid = fopen(lrcFile,'r');
line = fgetl(fid);
bin_size = str2double(line);
line = fgetl(fid);
fields = regexp(line,',','split');
data_chr_all = str2double(fields);
line = fgetl(fid);
fields = regexp(line,',','split');
data_bin_all = str2double(fields);
results = textscan(fid,repmat('%f',1,length(data_chr_all)),'Delimiter',',');
data_lrc_all = cell2mat(results);
clear results;
fclose(fid);

data_rc_all = 2.^data_lrc_all;
cluster_ids = unique(cell_labels);
for c = 1:length(cluster_ids)
    tv = cell_labels == cluster_ids(c);
    data_rc_all(tv,:) = data_rc_all(tv,:)/(median(reshape(data_rc_all(tv,:),[],1))+eps);
end
% for i = 1:size(data_rc_all,1)
%     data_rc_all(i,:) = data_rc_all(i,:)/(median(data_rc_all(i,:))+eps);
% end
data_lrc_all = log2(data_rc_all+eps);

chromosomes = reshape(unique(data_chr_all),1,[]);
sorted_indxs = [];
for i = 1:length(chromosomes)
    indxs = find(data_chr_all == chromosomes(i));
    sorted_indxs = [sorted_indxs indxs];
end

data_lrc_all = data_lrc_all(:,sorted_indxs);
data_chr_all = data_chr_all(sorted_indxs);
data_bin_all = data_bin_all(sorted_indxs);

end