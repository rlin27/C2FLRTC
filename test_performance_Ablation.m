clear all
close all
filename=cell(8,1);
filename{1}='airplane';
filename{2}='baboon';
filename{3}='barbara';
filename{4}='facade';
filename{5}='house';
filename{6}='lena';
filename{7}='peppers';
filename{8}='sailboat';
filename{9}='2';
filename{10}='5';
filename{11}='9';
filename{12}='17';
rng(602)
for ff=6:6
    myName=sprintf('TestImages/%s.bmp',filename{ff});
    A=imread(myName);
    A=imresize(A,[256,256]);
    figure(1);  imshow(A);
    
    myrate=0.90:0.1:0.90;
    %myrate=[myrate,0.99];
    myResult=cell(2,numel(myrate));
    A=double(A)/255.0;
    
    for iterate=1:numel(myrate)
        tic
        
        rate=1 - myrate(iterate);
        [row, col, channel]=size(A);
        B=zeros([row, col, channel]);
        mark=true([row, col, channel]);
        
        counter=1;
        for i=1:row
            for j=1:col
                for k=1:channel
                    if(rand()<rate)
                        index(counter,1)=i;
                        index(counter,2)=j;
                        index(counter,3)=k;
                        value(counter)=A(i,j,k);
                        B(i,j,k)=A(i,j,k);
                        mark(i,j,k)=false;
                        counter=counter+1;
                    end
                end
            end
        end
        figure(2);imshow(B);
        
        tsize=[row, col, channel];
        N=3;
        lambda=0.02;
        alpha=[1/N, 1/N, 1/N];
        beta=[1,1,0];
        
        fprintf('--------------TR_LRTV2-------------------\n');
        lambda_1=0.5;
        lambda_2=1000;%1000
        
        scale_num=4;
        overlap_pixel=5; %a half
        thre=0.4; %to decide whether replace the particular patch
        rse_r{ff,iterate}=[];
        r_n{ff,iterate}=[];
        psnr_r{ff,iterate}=[];
        diff_r{ff,iterate}=[];
        
        % padding A and B
        A_pad=zeros(row+2*overlap_pixel, col+2*overlap_pixel, channel);
        A_pad(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=A;
        B_pad=zeros(row+2*overlap_pixel, col+2*overlap_pixel, channel);
        B_pad(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=B;
        B_padimage{ff,iterate}=B_pad;
        % first complete the whole image
        [index_temp]=find(B_pad(:));
        value_temp=B_pad(index_temp);
        Z_TRLRTV2_temp=LRTC_TV_II(index_temp,value_temp, lambda_1, lambda_2 ,alpha, beta, size(B_pad), N ,300);
        Z_LRTCTV2_completewhole{ff,iterate}=Z_TRLRTV2_temp;
                save data Z_TRLRTV2_temp
%          load data
        tmpRSE_LRTC_TV2{ff,iterate}(1,1)=rse(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
        tmpPSNR_LRTC_TV2{ff,iterate}(1,1)=psnr(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
        rse_r{ff,iterate}=[rse_r{ff,iterate};tmpRSE_LRTC_TV2{ff,iterate}(1,1)];
        r_n{ff,iterate}=[r_n{ff,iterate};1];
        psnr_r{ff,iterate}=[psnr_r{ff,iterate}; tmpPSNR_LRTC_TV2{ff,iterate}(1,1)];
        diff_r{ff,iterate}=[diff_r{ff,iterate};0];
        
        
        
        figure(6);
        imshow(Z_TRLRTV2_temp);
        % do patch_wise completion
        id=1;
        for i=2:scale_num
            if i==2
                thre=thre;
            else
                thre=max(diff_r{ff,iterate}(1:5));
            end
            
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
%                     B_temp=Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                    extra_B=Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                    ra=randperm(prod(size(extra_B)));
                    part_ra=ra(1:floor(0.5*prod(size(extra_B))*(1-myrate(iterate))));
                    part_ra=setdiff(part_ra,index_temp);
                    B_temp(part_ra)=extra_B(part_ra).*(1+0.0*randn(size(part_ra)));
                    [index_temp]=find(B_temp(:));
                    value_temp=B_temp(index_temp);

                    if i<4
                        lambda_1_temp=lambda_1/((2^i));
                        lambda_2_temp=lambda_2/((2^i)^2); % 100 for 2*2;
%                         lambda_1_temp=lambda_1/((2*2));
%                         lambda_2_temp=lambda_2/((2*2)); % 100 for 2*2;
                    else % thre should be changed also, 0.4 would be ok
                        lambda_1_temp=lambda_1/((2^3));
                        lambda_2_temp=lambda_2/((2^(3))^2);
%                         lambda_1_temp=lambda_1/((2*2));
%                         lambda_2_temp=lambda_2/((2*(2)));
                        
                        B_temp_large=imresize(B_temp,[floor(size(B_temp,1)*2^(i-3)),floor(size(B_temp,2)*2^(i-3))],'box');%'box','nearest'
%                         B_temp_large=imresize(B_temp,[row/(2^(3-1))+2*overlap_pixel, col/(2^(3-1))+2*overlap_pixel],'box');
                        
                        [index_temp_large]=find(B_temp_large(:));
                        value_temp_large=B_temp_large(index_temp_large);
                        patch_size_temp=size(B_temp_large);
                    end
                    
                    if i<4
                        Z_TRLRTV2_patch=LRTC_TV_II(index_temp,value_temp, lambda_1_temp, lambda_2_temp ,alpha, beta, patch_size, N ,100);
                    else
                        Z_TRLRTV2_patch=LRTC_TV_II(index_temp_large,value_temp_large, lambda_1_temp, lambda_2_temp ,alpha, beta, patch_size_temp, N ,100);
                        Z_TRLRTV2_patch=imresize(Z_TRLRTV2_patch,[size(B_temp,1),size(B_temp,2)]);
                    end
                    
                    diff=rse(Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:),Z_TRLRTV2_patch);
                    diff_r{ff,iterate}=[diff_r{ff,iterate};diff];
                    % compare the patchs
                    figure(20);
                    subplot(3,1,1);
                    imshow(A_pad((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:))
                    subplot(3,1,2);
                    imshow(Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:))
                    subplot(3,1,3);
                    imshow(Z_TRLRTV2_patch);
                    
                    if diff<thre
                        temp=Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                        
                        Z_TRLRTV2_temp((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:)=Z_TRLRTV2_patch.*conv_mat1+temp.*conv_mat2;
                        r_n{ff,iterate}=[r_n{ff,iterate};1];
                    else
                        r_n{ff,iterate}=[r_n{ff,iterate};0];
                    end
                    figure(21);
                    imshow(Z_TRLRTV2_temp)
                    
%                     fprintf('rse with true patch\n')
                    tmpRSE_LRTC_TV2{ff,iterate}(1,i)=rse(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
%                     fprintf('psnr with true patch\n')
                    tmpPSNR_LRTC_TV2{ff,iterate}(1,i)=psnr(A,Z_TRLRTV2_temp(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:));
                    rse_r{ff,iterate}=[rse_r{ff,iterate};tmpRSE_LRTC_TV2{ff,iterate}(1,i)];
                    psnr_r{ff,iterate}=[psnr_r{ff,iterate};tmpPSNR_LRTC_TV2{ff,iterate}(1,i)];
                    
                    
                end
            end
%                                 figure(23);
%                     imshow(Z_TRLRTV2_temp)
            
        end
        Z_LRTCTV2_completepatch{ff,iterate}=Z_TRLRTV2_temp;
        
%         save simple_result_noTV rse_r r_n psnr_r diff_r Z_LRTCTV2_completepatch Z_LRTCTV2_completewhole B_padimage tmpRSE_LRTC_TV2 tmpPSNR_LRTC_TV2
        fprintf('per image per miss rate costs time:   ')
        toc
    end
end

% psnr=zeros(85,1);
% for i =1:8
%     psnr=psnr+psnr_r{i};
% end
% psnr=psnr/8;

%  for i=9:12
%  temp= Z_TRLRTV2_temp
%  figure(1);imshow(temp(6:end-5,6:end-5,:),'border','tight','initialmagnification','fit')
%  set (gcf,'Position',[0,0,500,500]);
%  axis normal;
%   end


% for i=1:8
%  temp= B;
%  figure(1);imshow(B,'border','tight','initialmagnification','fit')
%  set (gcf,'Position',[0,0,500,500]);
%  axis normal;
%  end