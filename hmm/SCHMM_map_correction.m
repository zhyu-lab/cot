function corrected_rc = SCHMM_map_correction(data_rc, data_map)
% Mappability correction for read counts

corrected_rc = data_rc;
m_all_map = median(data_rc,1);
int_map = floor(data_map*100);
int_map_u = unique(int_map);
m_map = zeros(length(int_map_u),size(data_rc,2));

for i = 1:length(int_map_u)
    tv = int_map == int_map_u(i);
    m_map(i,:) = median(data_rc(tv,:),1);
    corrected_rc(tv,:) = data_rc(tv,:).*repmat(m_all_map,sum(tv),1)./(repmat(m_map(i,:),sum(tv),1)+eps);
end

end