function corrected_rc = SCHMM_gc_correction(data_rc,data_gc)
% GC correction for read counts

corrected_rc = data_rc;
m_all_gc = median(data_rc,1);
int_gc = floor(data_gc*100);
int_gc_u = unique(int_gc);
m_gc = zeros(length(int_gc_u),size(data_rc,2));

for i = 1:length(int_gc_u)
    tv = int_gc == int_gc_u(i);
    m_gc(i,:) = median(data_rc(tv,:),1);
    corrected_rc(tv,:) = data_rc(tv,:).*repmat(m_all_gc,sum(tv),1)./(repmat(m_gc(i,:),sum(tv),1)+eps);
end

end