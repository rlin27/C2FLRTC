%---- Clear the Workspace ----
clear all
close all

%---- Download the test images ----
filename=cell(8,1);
filename{1}='airplane';
filename{2}='baboon';
filename{3}='barbara';
filename{4}='facade';
filename{5}='house';
filename{6}='lena';
filename{7}='peppers';
filename{8}='sailboat';

%---- Fix the random seed ----
rng(602)

%---- C2F ----
for ff=1:8  % load the test images
    myName=sprintf('test_images/%s.bmp',filename{ff});  % get the path of the selected image
    A=imread(myName);  % read the original image
    figure(1); imshow(A); title('Original Figure'); % print the orignal image
    
    myrate=0.70:0.1:0.90;  % set the missing rate
    myResult=cell(2,numel(myrate));  % build a cell to store the results
    A=double(A)/255.0;  % normalize the values of the image
    
    for iterate=1:numel(myrate)
        tic  % start the timer
        
        rate=1 - myrate(iterate);  % the observed rate
        [row, col, channel]=size(A);  % record the size of the image
        B=zeros([row, col, channel]);  % build a all-zero tensor having the same size as A
        mark=true([row, col, channel]);  % build a all-one tensor having the same size as A
        
        counter=1;  
        for i=1:row  % for each row 
            for j=1:col  % for each column
                for k=1:channel  % for each channel
                    if(rand()<rate)  % if the random number is smaller than the rate
                        index(counter,1)=i;  % each culumn of index records a coordinate of a pixel
                        index(counter,2)=j;
                        index(counter,3)=k;
                        value(counter)=A(i,j,k);  % each element of value records the value of the pixel
                        B(i,j,k)=A(i,j,k);  % B is the image with missing pixels
                        mark(i,j,k)=false;  % mark is a tensor with binary values
                        counter=counter+1;  % the number of observed pixels
                    end
                end
            end
        end
        figure(2);imshow(B); title(['Figure with missing ratio ', num2str(1-rate)]); % show the image with missing value
        
        tsize=[row, col, channel];  % the size of the image
        
        % hyper parameters
        tau = [32, 32, 1];
        param.incR{1} = [1, 2, 4, 8, 16, 24, 32];
        param.incR{2} = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 192, 256];
        param.incR{3} = [1, 2, 4, 8, 16, 24, 32];
        param.incR{4} = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 192, 256];
        param.incR{5} = [1];
        param.incR{6} = [1, 3];
        param.delta = 1e-5;
        sc = 255;
        
        fprintf('-------------- SPC -------------------\n');
        
        scale_num=5;
        overlap_pixel=5;  % a half
        thre=0.1;  % to decide to replace the particular patch or not
        
        % initialize the metrics
        rse_r{ff,iterate}=[];  % RSE
        r_n{ff,iterate}=[];  % replace vector
        psnr_r{ff,iterate}=[];  % PSNR
        diff_r{ff,iterate}=[];  % differences between the coarse image and the restored patch
        
        % padding A and B
        A_pad=zeros(row+2*overlap_pixel, col+2*overlap_pixel, channel);
        A_pad(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=A;
        B_pad=zeros(row+2*overlap_pixel, col+2*overlap_pixel, channel);
        B_pad(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=B;
        B_padimage{ff,iterate}=B_pad;
        % first complete the whole image
        [index_temp]=find(B_pad(:));
        value_temp=B_pad(index_temp);
        
        %%%%%%%%%% STDC completion %%%%%%%%%%
        Qms = (B_pad ~= 0);
        Z_TRLRTV2_temp=MDT_Tucker_incR(B_pad,Qms, tau, param);
        
        Z_LRTCTV2_completewhole{ff,iterate}=Z_TRLRTV2_temp;

        tmpRSE_LRTC_TV2{ff,iterate}(1,1)=rse(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
        tmpPSNR_LRTC_TV2{ff,iterate}(1,1)=psnr(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
        
        % SSIM
        tmpSSIM_LRTC_TV2{ff,iterate}(1,1)=calc_ssim(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
        
        rse_r{ff,iterate}=[rse_r{ff,iterate}; tmpRSE_LRTC_TV2{ff,iterate}(1,1)];
        r_n{ff,iterate}=[r_n{ff,iterate};1];
        psnr_r{ff,iterate}=[psnr_r{ff,iterate}; tmpPSNR_LRTC_TV2{ff,iterate}(1,1)];
        diff_r{ff,iterate}=[diff_r{ff,iterate};0];
        
        % save the result got by pure LRTC
        store_name=sprintf('pure_SPC_restore_%s.mat',filename{ff});
        save(store_name,'Z_TRLRTV2_temp')
        
        % coarse level completion
        figure(4);
        imshow(Z_TRLRTV2_temp); title('Image Completed by pure LRTC method')
        
        % do patch-wise completion
        for i=2:scale_num
            patch_size=[row/(2^(i-1))+2*overlap_pixel, col/(2^(i-1))+2*overlap_pixel, channel];
            conv_mat1=1/2*ones(patch_size);
            conv_mat1(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=ones(patch_size(1)-2*overlap_pixel,patch_size(2)-2*overlap_pixel,patch_size(3));
            conv_mat2=1/2*ones(patch_size);
            conv_mat2(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=zeros(patch_size(1)-2*overlap_pixel,patch_size(2)-2*overlap_pixel,patch_size(3));
            
            for row_wise=1:2^(i-1)
                for col_wise=1:2^(i-1)
                    B_temp=B_pad((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                    [index_temp]=find(B_temp(:));
                    value_temp=B_temp(index_temp);
                     B_temp=Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                    sizeBtemp=prod(size(B_temp));
                    vecBtemp=randperm(sizeBtemp);
                    index_temp=[index_temp;vecBtemp(1:floor(length(vecBtemp)*0.1))'];
                    value_temp=B_temp(index_temp);
                    
                    if i<6
                        %%%%%%%%%% STDC completion %%%%%%%%%%
                        Qms = (B_pad ~= 0);
                        Z_TRLRTV2_temp=MDT_Tucker_incR(B_pad,Qms, tau, param);
                    else
                        %%%%%%%%%% STDC completion %%%%%%%%%%
                        Qms = (B_pad ~= 0);
                        Z_TRLRTV2_temp=MDT_Tucker_incR(B_pad,Qms, tau, param);
                        
                        Z_TRLRTV2_patch=imresize(Z_TRLRTV2_patch,[size(B_temp,1),size(B_temp,2)]);
                    end
                    
                    diff=rse(Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:),Z_TRLRTV2_patch);
                    diff_r{ff,iterate}=[diff_r{ff,iterate};diff];
                    
                    % compare the patchs
                    figure(5);
                    subplot(3,1,1);
                    imshow(A_pad((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:))
                    title('Original Patch');
                    subplot(3,1,2);
                    imshow(Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:))
                    title('k-th Completed Patch');
                    subplot(3,1,3);
                    imshow(Z_TRLRTV2_patch);
                    title('(k+1)-th Completed Patch'); 
                    
                    if diff<thre
                        temp=Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                        
                        Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:)=Z_TRLRTV2_patch.*conv_mat1+temp.*conv_mat2;
                        r_n{ff,iterate}=[r_n{ff,iterate};1];
                    else
                        r_n{ff,iterate}=[r_n{ff,iterate};0];
                    end
                    figure(6);
                    imshow(Z_TRLRTV2_temp); title('Temporary Completed Image');
                    
                    tmpRSE_LRTC_TV2{ff,iterate}(1,i)=rse(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
                    tmpPSNR_LRTC_TV2{ff,iterate}(1,i)=psnr(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
                    
                    % SSIM
                    tmpSSIM_LRTC_TV2{ff,iterate}(1,i)=calc_ssim(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
                    
                    rse_r{ff,iterate}=[rse_r{ff,iterate};tmpRSE_LRTC_TV2{ff,iterate}(1,i)];
                    psnr_r{ff,iterate}=[psnr_r{ff,iterate};tmpPSNR_LRTC_TV2{ff,iterate}(1,i)];
                end
            end
        end
        Z_LRTCTV2_completepatch{ff,iterate}=Z_TRLRTV2_temp;
        
        save C2F_SPC_Results rse_r r_n psnr_r diff_r Z_LRTCTV2_completepatch Z_LRTCTV2_completewhole B_padimage tmpRSE_LRTC_TV2 tmpPSNR_LRTC_TV2
        fprintf('per image per miss rate costs time:   ')
        toc
    end
end