%%
% for grappa_v2, data dimension order is [PE, RO, CHA]
%%
clear;

load('brain.mat')
kspace2D = ksp;

%% Generate downsampling and estimate CSM
rfac = 4;
[d1,d2,d3] = size(kspace2D); 
ndim = d1;    %phase encoding direction
off = 0;       %starting sampling location

nencode = 36;  % The number of ACS lines
 
num_block = 3;
num_column = 5;

acs_line_loc = (ndim/2+1-nencode/2):(ndim/2+nencode/2);
pe_loc = (off+1):rfac:(d1-off);

acq_idx = zeros(d1,1);
acq_idx(pe_loc) = 1;
acq_idx(acs_line_loc) = 1;
NetR = d1 / sum(acq_idx)

%% GRAPPA recon
sp = zeros(d1,d2);
sp(acs_line_loc,:) = 2;
csm = ismrm_estimate_csm(kspace2D,sp);

k_space_red = kspace2D(pe_loc,:,:);
acs_data = kspace2D(acs_line_loc,:,:);
tic
[full_fourier_data0] = grappa_v2(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column);
toc
if size(full_fourier_data0,1) < size(kspace2D,1)
    kspace2D_recon = zeros(size(kspace2D));
    kspace2D_recon(1:1:size(full_fourier_data0,1),:,:) = full_fourier_data0;
end
if size(full_fourier_data0,1) > size(kspace2D,1)
    kspace2D_recon = full_fourier_data0(1:1:size(kspace2D,1),:,:);
end
if size(full_fourier_data0,1) == size(kspace2D,1);
    kspace2D_recon = full_fourier_data0;
end
im_recon = sum(ismrm_transform_kspace_to_image(kspace2D_recon,[1,2]).*conj(csm),3);
as(im_recon)

%% Direct g-factor comparison between GRAPPA and VCC-GRAPPA
kspace2D_tmp = permute(kspace2D,[3,2,1]);
kspace2D_vcc = VCC_signal_creation(kspace2D_tmp);
kspace2D_vcc = permute(kspace2D_vcc,[3,2,1]);

csm_vcc = ismrm_estimate_csm(kspace2D_vcc,sp);

%% VCC-grappa recon
k_space_red = kspace2D_vcc(pe_loc,:,:);
acs_data = kspace2D_vcc(acs_line_loc,:,:);
 
tic
[full_fourier_data0, ImgRecon0, coef0] = grappa_v2(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column);
% times_comp = 2;
% [full_fourier_data0, ImgRecon0, coef0] = nonlinear_grappa(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column,times_comp);
toc
if size(full_fourier_data0,1) < size(kspace2D,1)
    kspace2D_vcc_recon = zeros(size(kspace2D_vcc));
    kspace2D_vcc_recon(1:1:size(full_fourier_data0,1),:,:) = full_fourier_data0;
end
if size(full_fourier_data0,1) > size(kspace2D,1)
    kspace2D_vcc_recon = full_fourier_data0(1:1:size(kspace2D,1),:,:);
end
if size(full_fourier_data0,1) == size(kspace2D,1);
    kspace2D_vcc_recon = full_fourier_data0;
end
im_recon_vcc = sum(ismrm_transform_kspace_to_image(kspace2D_vcc_recon,[1,2]).*conj(csm_vcc),3);
as(im_recon_vcc)


im_true_coil = ismrm_transform_kspace_to_image(kspace2D,[1,2]);
im_true = sum(im_true_coil .* conj(csm),3);

im_true_vcccoil = ismrm_transform_kspace_to_image(kspace2D_vcc,[1,2]);
im_true_vcc = sum(im_true_vcccoil .* conj(csm_vcc),3);

im_diff_grappa = mat2gray(abs(im_true)) - mat2gray(abs(im_recon));
rmse_grappa = norm(im_diff_grappa(:))/norm(im_true(:))

im_diff_vccgrappa = mat2gray(abs(im_true)) - mat2gray(abs(im_recon_vcc));
rmse_vccgrappa = norm(im_diff_vccgrappa(:))/norm(im_true(:))


is_pseudo_replica = 0;
if (is_pseudo_replica)
    noise_level = 0.001*max(abs(im_recon_vcc(:)));
    
    reps = 500;
    im_size = [size(kspace2D,1),size(kspace2D,2)];
    img_noise_rep = zeros([im_size,reps]);
    img_vcc_noise_rep = zeros([im_size,reps]);
    noise_rep = zeros([im_size,reps]);
    ref_img_noise_rep = zeros([im_size,reps]);
    tic
    parfor r = 1:reps
        % white gaussian noise scaled by noise_level
        noise_white = noise_level*complex(randn(size(kspace2D)),randn(size(kspace2D)));
        
        % add noise to image
        im_ref_noise = im_true_coil + noise_white;
        
        % simulate noisy k-space data
        data_noise = ismrm_transform_image_to_kspace(im_ref_noise,[1,2]);
        
        % simulate VCC-GRAPPA downsampling
        kspace2D_tmp = permute(data_noise,[3,2,1]);
        kspace2D_vcc = VCC_signal_creation(kspace2D_tmp);
        kspace2D_vcc = permute(kspace2D_vcc,[3,2,1]);
        data_noise = kspace2D_vcc;
        
        im_ref_vcc_noise = ismrm_transform_kspace_to_image(kspace2D_vcc,[1,2]);
        
        
        k_space_red = data_noise(pe_loc,:,:);
        % acs_data = data_noise(acs_line_loc,:,:);
        
        % GRAPPA reconstruction
%                 tic
%                 [full_fourier_data0, ImgRecon0, coef0] = grappa_v2(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column);
%                 toc
        % nonlinear grappa reconstruction
        tic
        times_comp = 2;
        [full_fourier_data0, ImgRecon0, coef0] = nonlinear_grappa(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column,times_comp);
        toc
%         
        if size(full_fourier_data0,1) < size(kspace2D,1)
            kspace2D_vcc_recon1 = zeros(size(kspace2D));
            kspace2D_vcc_recon1(1:1:size(full_fourier_data0,1),:,:) = full_fourier_data0;
            
        end
        if size(full_fourier_data0,1) > size(kspace2D,1)
            kspace2D_vcc_recon1 = full_fourier_data0(1:1:size(kspace2D,1),:,:);
        end
        if size(full_fourier_data0,1) == size(kspace2D,1);
            kspace2D_vcc_recon1 = full_fourier_data0;
        end
        im_recon = sum(ismrm_transform_kspace_to_image(kspace2D_vcc_recon1,[1,2]).*conj(csm_vcc),3);
        
        % recorded image
        img_noise_rep(:,:,r) = im_recon;
        ref_img_noise_rep(:,:,r) = sum(im_ref_noise.*conj(csm),3);
        img_vcc_noise_rep(:,:,r) = sum(im_ref_vcc_noise.*conj(csm_vcc),3);
        noise_rep(:,:,r) = sum(noise_white.*conj(csm),3);
    end
    toc
    
    img_noise_rep = reshape(img_noise_rep,[im_size,reps]);
    rep_dim = length(size(img_noise_rep));
    
    std_pseudo = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
    % std_pseudo = std(abs(img_noise_rep),[],rep_dim);))))
    
    std_noise = std(abs(noise_rep + max(abs(noise_rep(:)))),[],rep_dim);
    
    % std_full = std(abs(ref_img_noise_rep + max(abs(ref_img_noise_rep(:)))),[],rep_dim);
    % std_full = std(abs(ref_img_noise_rep),[],rep_dim);
    std_full = std(abs(img_vcc_noise_rep),[],rep_dim);
    
    gmap_pseudo = std_pseudo ./(std_full.*sqrt(rfac));    % Be careful about the scaling factor sqrt(NetR)
    
    gmap_pseudo(gmap_pseudo < eps) = 1;
    
    as(gmap_pseudo)
    sprintf('VCC GRAPPA g_mean by Pseudo Replica : %f, g_max : %f',mean(gmap_pseudo(:)),max(gmap_pseudo(:)))
end