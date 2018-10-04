function bdata = baseline_correction(data,indices)
% data is of dim chans x time
% indices refer to the second dim
bdata     = data - repmat(mean(data(:,indices),2), [1,size(data,2)]);