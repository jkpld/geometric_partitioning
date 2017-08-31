function locs = max_curvature_idx(kappa)

if iscell(kappa)
    locs = cell(1,numel(kappa));
    for i = 1:numel(kappa)
        locs{i} = curvature_max__(kappa{i});
    end
else
    locs = curvature_max__(kappa);
end

end


function locs = curvature_max__(kappa)
    
    MAD_multiplier = 2;
    threshold = max(0, MAD_multiplier*mad(kappa,1) + median(kappa,'omitnan'));
    
    H = 3;
    
    nan_idx = isnan(kappa);
    cntrNum = cumsum(nan_idx)+1;
    cntrNum(nan_idx) = 0;
    idx = (1:numel(kappa))';
    
    locs = cell(1,cntrNum(end));
    
    for i = 1:cntrNum(end)
        inds = cntrNum == i;
        c_idx = idx(inds);
        c_kappa = kappa(inds);

        c_idx = padarray(c_idx,floor(H/2)+1, 'circular');
        c_kappa = padarray(c_kappa,floor(H/2)+1, 'circular');
        
        c_locs = (c_kappa == imdilate(c_kappa,ones(H,1))) & (c_kappa > threshold);
        locs{i} = unique(c_idx(c_locs));
        
    end
    
    locs = cat(1,locs{:});
end
