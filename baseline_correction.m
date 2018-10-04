function bdata = baseline_correction(rdata,indices)
% data is of dim chans x time
% indices refer to the second dim
if isempty(indices)
    warning('W: baseline correction needs a window defined by either its limits or by all members.');
elseif length(indices)==2 
    indices = indices(1):indices(2);
end
bdata     = rdata - repmat(mean(rdata(:,indices),2), [1,size(rdata,2)]);