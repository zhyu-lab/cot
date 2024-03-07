function plot_cn_results(result_dir)

label_file = [result_dir '/label.txt'];
fid = fopen(label_file,'r');
line = fgetl(fid);
barcodes = regexp(line,',','split');
line = fgetl(fid);
cell_labels = str2double(regexp(line,',','split'));
cluster_ids = unique(cell_labels);
fclose(fid);

lrc_file = [result_dir '/lrc.txt'];
[data_lrc_all,data_chr_all,data_bin_all,bin_size] = load_data(lrc_file);
data_spos_all = (data_bin_all-1)*bin_size+1;
data_epos_all = data_bin_all*bin_size;

fid = fopen([result_dir '/paras.csv'],'r');
if fid == -1
    error(['Can not open file: ' result_dir '/paras.csv']);
end
results = textscan(fid,'%f%*f','headerlines',1,'delimiter',',');
fclose(fid);
o_all = results{1};

fid = fopen([result_dir '/ploidy.csv'],'r');
if fid == -1
    error(['Can not open file: ' result_dir '/ploidy.csv']);
end
line = fgetl(fid);
fields = regexp(line,',','split');
acn_all = str2double(fields);
fclose(fid);

fid = fopen([result_dir '/segments.csv'],'r');
if fid == -1
    error(['Can not open file: ' result_dir '/segments.csv']);
end
results = textscan(fid,'%f%f%f%f%f','headerlines',1,'delimiter',',');
fclose(fid);
cell_seg_all = results{1};
chr_seg_all = results{2};
spos_seg_all = results{3};
epos_seg_all = results{4};
cn_seg_all = results{5};

fid = fopen([result_dir '/copynumber.csv'],'r');
if fid == -1
    error(['Can not open file: ' result_dir '/copynumber.csv']);
end
line = fgetl(fid);
bin_count = length(regexp(line,',','split'));
results = textscan(fid,repmat('%f',1,bin_count),'delimiter',',');
fclose(fid);
cn_bin_all = cell2mat(results);

line_style = 'k-';
LineWidth = 0.5;
MarkerSize = 4;
FontSize = 12;

lcr_colors = [0.9 0.9 0.9;
    0 0.9 0;
    0 0 0.9;
    0.9 0 0];

chromosomes = unique(data_chr_all);
max_pos = zeros(1,length(chromosomes));
for c = 1:length(chromosomes)
    tv = data_chr_all == chromosomes(c);
    max_pos(c) = max(data_epos_all(tv));
end
ratio = max_pos/sum(max_pos);
xtick = cumsum([0 ratio(1:end-1)])+ratio/2;

figure(1);
clf;

u_labels = unique(cell_labels);
rows = 4;
cols = ceil(length(u_labels)/(rows/2));

for k = 1:length(u_labels)
    disp(num2str(sum(cell_labels == u_labels(k))))
    cell_indxs = find(cell_labels == u_labels(k));
    tmp = [];
    for i = 1:length(cell_indxs)
        cn_seg = cn_bin_all(cell_indxs(i),:);
        tmp = [tmp; cn_seg];
    end
    m_cns = mean(tmp,1);
    dists = sum((tmp-repmat(m_cns,length(cell_indxs),1)).^2,2);
    [~, j] = min(dists);
    cell_id = cell_indxs(j);
    
    tv = cell_seg_all == cell_id;
    chr_seg = chr_seg_all(tv);
    pstart_seg = spos_seg_all(tv);
    pend_seg = epos_seg_all(tv);
    cn_seg = cn_seg_all(tv);
    acn = acn_all(cell_id);
    
    tv = cell_labels(cell_id) == cluster_ids;
    o = o_all(tv);

    data_lrc_cell = data_lrc_all(cell_id,:);
    subplot(rows,cols,floor((k-1)/cols)*cols*2+rem(k-1,cols)+1);
    set(gca,'YGrid','on','Box','on','FontSize',FontSize);
    hold on
    pre_x = 0;
    chr_epos = zeros(length(chromosomes),1);
    for c = 1:length(chromosomes)
        tv = data_chr_all == chromosomes(c);
        data_lcr = data_lrc_cell(tv);
        data_spos = data_spos_all(tv);
        data_epos = data_epos_all(tv);
        x = data_epos*ratio(c)/max_pos(c)+pre_x;
        indx1 = find(chr_seg == chromosomes(c));
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            tv = data_spos >= pstart_seg(j) & data_epos <= pend_seg(j);
            if sum(tv) == 0
                continue;
            end
            if CN < 1
                m = 1;
            else
                m = CN+1;
            end
            if m > 4
                m = 4;
            end
            plot(x(tv),data_lcr(tv),'.','MarkerSize',MarkerSize, 'Color', lcr_colors(m,:));
        end
        % plot expected LCR mean values
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            if CN == 0
                CN = 0.001;
            end
            lrc_mean = log2(CN/2)+o;
            indx = find(data_spos >= pstart_seg(j) & data_epos <= pend_seg(j));
            if isempty(indx)
                continue;
            end
            plot([x(indx(1)) x(indx(end))],[lrc_mean lrc_mean],'k-','LineWidth',1.5);
        end
        chr_epos(c) = max(x);
        pre_x = pre_x+ratio(c);  
    end
    for c = 1:length(chromosomes)-1
        plot([chr_epos(c) chr_epos(c)],[-3 3],line_style,'LineWidth',LineWidth)
    end
    ylabel('LRC');
    set(gca,'ytick',[-2 0 2])
    axis([0 1 -2.5 2.5])
%     tmp = [barcode ', ploidy = ' num2str(acn)];
    tmp = ['Cluster ' num2str(k) ', ploidy = ' num2str(acn)];
	title(tmp,'fontsize',FontSize+2);
    set(gca,'XTick',[])
    
    subplot(rows,cols,floor((k-1)/cols)*cols*2+rem(k-1,cols)+cols+1);
    set(gca,'YGrid','on','Box','on','FontSize',FontSize);
    set(gca,'YGrid','on','Box','on','FontSize',FontSize);
%     set(gca,'YTick',[1:1:7],'Box','on');
%     set(gca,'YTickLabel',{'1','2','3','4','5','6','>=7'});
    set(gca,'YTick',[1:2:7],'Box','on');
    set(gca,'YTickLabel',{'1','3','5','>=7'});
    hold on
    pre_x = 0;
    for c = 1:length(chromosomes)
        tv = data_chr_all == chromosomes(c);
        data_lcr = data_lrc_cell(tv);
        data_spos = data_spos_all(tv);
        data_epos = data_epos_all(tv);
        x = data_epos*ratio(c)/max_pos(c)+pre_x;
        indx1 = find(chr_seg == chromosomes(c));
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            indx = find(data_spos >= pstart_seg(j) & data_epos <= pend_seg(j));
            if sum(tv) == 0
                continue;
            end
            if CN > 7
                CN = 7;
            end
            plot([x(indx(1)) x(indx(end))],[CN CN],'-r','LineWidth',2.0);
        end
        pre_x = pre_x+ratio(c); 
    end
    for c = 1:length(chromosomes)-1
        plot([chr_epos(c) chr_epos(c)],[-0.1 7.5],line_style,'LineWidth',LineWidth)
    end
    ylabel('Copy number');
    xlabel('Chromosome');
    axis([0 1 0.5 7.5])
    set(gca,'XTick',xtick(1:4:end));
    set(gca,'XTickLabel',mat2cell(chromosomes(1:4:end),1,length(chromosomes(1:4:end))));
    
end

%save figure
figpath = [result_dir '/segment.png'];
eval(['print -dpng -r300 ' figpath])

end

function [data_lrc_all,data_chr_all,data_bin_all,bin_size] = load_data(lrc_file)

fid = fopen(lrc_file, 'r');
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
for i = 1:size(data_rc_all,1)
    data_rc_all(i,:) = data_rc_all(i,:)/(median(data_rc_all(i,:))+eps);
end
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

