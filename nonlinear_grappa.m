function [full_fourier_data, rec_img, coef0] = nonlinear_grappa(reduced_fourier_data, ORF, pe_loc, acs_data, acs_line_loc, num_block, num_column, times_comp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Yuchou Chang, University of Wisconsin - Milwaukee
% Email: yuchou@uwm.edu; leiying@uwm.edu
% Citation: Y. Chang, D. Liang, L. Ying, "Nonlinear GRAPPA: A Kernal Approach 
%           to Parallel MRI Reconstruction". Magn. Reson. Med. 2012
% Created on Oct. 12, 2011

% Input parameters

% reduced_fourier_data : undersampled k-space data
% ORF : outer reduction factor
% pe_loc : undersampled phase-encoding lines' location
% acs_data : auto-calibration signal data (middle region of k-space)
% acs_line_loc : auto-calibration signal lines' location
% num_block : number of blocks
% num_column : number of columns
% times_comp : times of the number of the first-order terms (the number of 
% the second-order terms = time_comp X the number of the first-order terms)

% Output parameters

% full_fourier_data : reconstructed k-space (with ACS replacement)
% rec_img : reconstructed image
% coef0 : coefficients for reconstruction

% time_comp parameter
% As the parameter time_comp increases, relevant second-order terms are
% added for reconstruction. 
% When time_comp = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get dimensions and initialization
[d1_reduced,d2,num_coil] = size(reduced_fourier_data);
d1 = d1_reduced*ORF;
if ORF==3
    d1=d1_reduced*ORF-2;
end
if ORF==5
    d1=d1_reduced*ORF-4;
end
if ORF==6
    d1=d1_reduced*ORF-2;
end


%Decide which lines are possible lines to fit
all_acquired_line_loc = unique(sort([pe_loc, acs_line_loc]));
combined_fourier_data = zeros(d1,d2,num_coil);
combined_fourier_data(pe_loc,:,:) = reduced_fourier_data;
combined_fourier_data(acs_line_loc,:,:) = acs_data;
ind_first = find(all_acquired_line_loc == acs_line_loc(1));
ind_last = find(all_acquired_line_loc == acs_line_loc(end));

%Form the structure that indicates where lines are fitted
line_group = cell(num_block,ORF-1);
for s = ind_first:ind_last
    for mode = 1:num_block
        for offset = 1:ORF-1
            tentative_line_ind = [all_acquired_line_loc(s)-offset-(mode-1)*ORF : ORF : all_acquired_line_loc(s)-offset+(num_block-1)*ORF-(mode-1)*ORF];
            valid_flag = 1;
            for t = 1:num_block
                if isempty(find(all_acquired_line_loc == tentative_line_ind(t)))
                    valid_flag = 0;
                    break;
                end
            end
            if valid_flag == 1
                line_group{mode,offset} = unique([line_group{mode,offset}; [all_acquired_line_loc(s),tentative_line_ind] ], 'rows');
            end
        end
    end
end

%Solve for the weighting coefficients
fit_coef = zeros(num_coil,ORF-1,num_block,(times_comp+1)*num_block*num_coil*num_column);

for jj = 1:num_coil
    for offset = 1:ORF-1
      for mode = 1:num_block
          fit_mat = zeros(num_block*num_coil*num_column, d2*size(line_group{mode,offset},1) );
          target_vec = zeros(1,d2*size(line_group{mode,offset},1) );
         
          for nn = 1:size(line_group{mode,offset},1)
              temp_data = combined_fourier_data( line_group{mode,offset}(nn,[2:end]), :,:);
              temp_data = permute(temp_data,[1 3 2]);
              temp_data = reshape(temp_data, [num_block*num_coil,d2]);
              fit_mat((num_block*num_coil*floor(num_column/2)+1):(num_block*num_coil*ceil(num_column/2)), [1+(nn-1)*d2 : nn*d2]) = temp_data;
              target_vec([1+(nn-1)*d2 : nn*d2]) = combined_fourier_data( line_group{mode,offset}(nn,1), :,jj);
          end

          transfer_matrix=fit_mat((num_block*num_coil*floor(num_column/2)+1):num_block*num_coil*(floor(num_column/2)+1),:);
          column_label = [[floor(num_column/2):-1:1],[1:floor(num_column/2)]];
          for column_idx = 1:num_column-1
              if column_idx <= floor(num_column/2)
                  fit_mat(num_block*num_coil*(column_idx-1)+1:num_block*num_coil*column_idx,:) = ...
                      [transfer_matrix(:,column_label(column_idx)+1:end) transfer_matrix(:,1:column_label(column_idx))];
              else
                  fit_mat(num_block*num_coil*column_idx+1:num_block*num_coil*(column_idx+1),:) = ...
                      [transfer_matrix(:,end-column_label(column_idx)+1:end) transfer_matrix(:,1:end-column_label(column_idx))];
              end
          end
          
          fit_mat_dim = reshape(fit_mat, [num_block*num_coil,num_column,d2*size(line_group{mode,offset},1)]);
          fit_mat_dim = reshape(fit_mat, [num_block,num_coil,num_column,d2*size(line_group{mode,offset},1)]);
          
          % nonlinear transformation
          new_fit_mat = zeros((times_comp+1)*size(fit_mat,1),size(fit_mat,2));
          new_fit_mat(1:num_block*num_coil*num_column,:) = fit_mat;
                  
          idx_comp = 0;
          for idx_adj_1 = 0 : ceil(num_block/2)-1
          	  for idx_adj_2 = 0 : ceil(num_coil/2)-1
              	  for idx_adj_3 = 0 : ceil(num_column/2)-1
                  	  fit_mat_shift = circshift(fit_mat_dim,[idx_adj_1 idx_adj_2 idx_adj_3 0]);
                      fit_mat_shift = reshape(fit_mat_shift, [num_block,num_coil*num_column,d2*size(line_group{mode,offset},1)]);
                      fit_mat_shift = reshape(fit_mat_shift, [num_block*num_coil*num_column,d2*size(line_group{mode,offset},1)]);

                      new_fit_mat(num_block*num_coil*num_column*(idx_comp+1)+1:num_block*num_coil*num_column*(idx_comp+1)+num_block*num_coil*num_column,:) = fit_mat.*fit_mat_shift;
                      idx_comp = idx_comp+1;

                      if idx_comp >= times_comp
                      	  break;
                      end
                  end
                  if idx_comp >= times_comp
                  	  break;
                  end
              end
              if idx_comp >= times_comp
              	  break;
              end
          end

          %fit_coef(jj,offset,mode,:) = target_vec/new_fit_mat;
          fit_coef(jj,offset,mode,:) = inv((new_fit_mat.')'*new_fit_mat.')*(new_fit_mat.')'*target_vec.';
      end
    end
end
clear temp_data;

%Generate the missing lines using superpositions
candidate_fourier_data = zeros(d1,d2,num_coil,num_block);
for mode = 1:num_block
    candidate_fourier_data(:,:,:,mode) = combined_fourier_data;
end
for ss = 1:d1
    if isempty(find(pe_loc == ss))
        offset = mod(ss-1,ORF);
        for mode = 1:num_block
            tentative_line_ind = [ORF*floor((ss-1)/ORF)+1-(mode-1)*ORF : ORF : ORF*floor((ss-1)/ORF)+1+(num_block-1)*ORF-(mode-1)*ORF];
            if max(tentative_line_ind) <= d1 & min(tentative_line_ind) >= 1
                temp_data = combined_fourier_data(tentative_line_ind,:,:);
                temp_data = permute(temp_data,[1 3 2]);
                fit_mat=zeros(num_block*num_coil*num_column,d2);
                temp_data=reshape(temp_data,[num_block*num_coil,d2]);
                
                fit_mat((num_block*num_coil*floor(num_column/2)+1):num_block*num_coil*(floor(num_column/2)+1),:) = temp_data;
                column_label = [[floor(num_column/2):-1:1],[1:floor(num_column/2)]];
                for column_idx = 1:num_column-1
                    if column_idx <= floor(num_column/2)
                        fit_mat(num_block*num_coil*(column_idx-1)+1:num_block*num_coil*column_idx,:) = ...
                            [temp_data(:,column_label(column_idx)+1:end) temp_data(:,1:column_label(column_idx))];
                    else
                        fit_mat(num_block*num_coil*column_idx+1:num_block*num_coil*(column_idx+1),:) = ...
                            [temp_data(:,end-column_label(column_idx)+1:end) temp_data(:,1:end-column_label(column_idx))];
                    end
                end
                
                fit_mat_dim = reshape(fit_mat, [num_block*num_coil,num_column,d2]);
                fit_mat_dim = reshape(fit_mat, [num_block,num_coil,num_column,d2]);
                
                % nonlinear transformation
                new_fit_mat = zeros(size(fit_mat,1)*(times_comp+1), d2);
                new_fit_mat(1:num_block*num_coil*num_column,:) = fit_mat;
                
                idx_comp = 0;
                for idx_adj_1 = 0 : ceil(num_block/2)-1
                    for idx_adj_2 = 0 : ceil(num_coil/2)-1
                        for idx_adj_3 = 0 : ceil(num_column/2)-1
                            fit_mat_shift = circshift(fit_mat_dim,[idx_adj_1 idx_adj_2 idx_adj_3 0]);
                            fit_mat_shift = reshape(fit_mat_shift, [num_block,num_coil*num_column,d2]);
                            fit_mat_shift = reshape(fit_mat_shift, [num_block*num_coil*num_column,d2]);

                            new_fit_mat(num_block*num_coil*num_column*(idx_comp+1)+1:num_block*num_coil*num_column*(idx_comp+1)+num_block*num_coil*num_column,:) = fit_mat.*fit_mat_shift;
                            idx_comp = idx_comp+1;

                            if idx_comp >= times_comp
                                break;
                            end
                        end
                        if idx_comp >= times_comp
                            break;
                        end
                    end
                    if idx_comp >= times_comp
                        break;
                    end
                end
                
                for jj = 1:num_coil
                    candidate_fourier_data(ss,:,jj,mode) = (squeeze(fit_coef(jj,offset,mode,:))).'*new_fit_mat;
                end
            else
                candidate_fourier_data(ss,:,:,mode) = 0;
            end
        end
    end
end

%Use ACS lines to obtain the goodness-of-fit coefficients
gof_coef = zeros(num_coil,ORF-1,num_block);
for jj = 1:num_coil
   for offset = 1:ORF-1
       fit_mat =[];
       target_vec = [];
       for ss = 1:length(acs_line_loc)
           if mod(acs_line_loc(ss)-1,ORF) == offset
               valid_flag = 1;
               for mode = 1:num_block
                   if isempty(find(line_group{mode,offset}(:,1) == acs_line_loc(ss)))
                       valid_flag = 0;
                       break;
                   end
               end
               if valid_flag == 1
                   temp_mat = [];
                   for mode = 1:num_block
                       temp_mat = [temp_mat; candidate_fourier_data(acs_line_loc(ss),:,jj,mode)];
                   end
                   fit_mat = [fit_mat,temp_mat];
                   target_vec = [target_vec,combined_fourier_data(acs_line_loc(ss),:,jj)];
               end
           end
       end
       gof_coef(jj,offset,:) = target_vec/fit_mat;
   end
end

%Combine the data from different modes using goodness-of-fit
full_fourier_data = combined_fourier_data;
for ss = 1:d1
    if isempty(find(all_acquired_line_loc == ss))
        offset = mod(ss-1,ORF);
        for jj = 1:num_coil
            for mode = 1:num_block
                full_fourier_data(ss,:,jj) = full_fourier_data(ss,:,jj)+gof_coef(jj,offset,mode)*candidate_fourier_data(ss,:,jj,mode);
            end
        end
    end
end

% full_fourier_data(acs_line_loc,:,:) = acs_data;

%Image reconstruction using IFFT2 and sum-of-squares
if nargout > 1
    coil_img = fftshift(fftshift(ifft2(ifftshift(ifftshift(full_fourier_data,1),2)),1),2);
    rec_img = sqrt(sum(abs(coil_img).^2,3));
end

coef0=fit_coef;