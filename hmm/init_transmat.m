function T = init_transmat(num_block,block_size,clamp_thres)
%block diagonal transition matrix
%clamp_thres is used to make sure that
%the diagonal elements in each block are no less than it.
% transitions between blocks are not allowed.

num_elements = num_block*block_size;
T = zeros(num_elements,num_elements);
for k = 1:num_block
    indxs = (k-1)*block_size+1:k*block_size;
    T(indxs,indxs) = 1;
    for j = 1:length(indxs)
        i = indxs(j);
        temp = T(i,:);
        tmp = sum(temp);
        if T(i,i) < clamp_thres*tmp
            temp(i) = 0;
            T(i,:) = (1-clamp_thres).*temp/sum(temp);
            T(i,i) = clamp_thres;
        else
            T(i,:) = temp/(tmp+eps);%to avoid tmp = 0
        end
    end
end

end



