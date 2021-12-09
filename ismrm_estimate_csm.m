function csm = ismrm_estimate_csm(kspace2D,sp)
fprintf('Estimating coil sensitivity ... \n');
nc = size(kspace2D,3);
source_data = kspace2D .* repmat(sp > 1,[1,1,nc]);
data_mask = sp > 1;
f = hamming(max(sum(data_mask,1))) * hamming(max(sum(data_mask,2)))';
fmask = zeros(size(source_data));
fmask((1:size(f,1))+bitshift(size(source_data,1),-1)-bitshift(size(f,1),-1), ...
    (1:size(f,2))+bitshift(size(source_data,2),-1)-bitshift(size(f,2),-1), :) = ...
    repmat(f, [1 1 size(source_data,3)]);
csm = ismrm_transform_kspace_to_image(source_data .* fmask, [1 2]);
csm = ismrm_estimate_csm_walsh(csm); %Estimate coil sensitivity maps.

fprintf('done.\n');