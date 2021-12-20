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
ss

avg_rr_record{8, 4} = []
%---- C2F ----
for ff=1:8  % load the test images
    myName=sprintf('test_images/%s.bmp',filename{ff});  % get the path of the selected image
    A=imread(myName);  % read the original image
    figure(1); imshow(A); title('Original Figure'); % print the orignal image
    
    myrate=0.90;  % set the missing rate
    myResult=cell(2,numel(myrate));  % build a cell to store the results
    A=double(A)/255.0;  % normalize the values of the image
    
    A_rr = relative_rank(A);
    avg_rr_record{ff,1} = A_rr / 255;
    
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
        
        tsize=[row, col, channel];  % the size of the image
        
        % hyper parameters
        scale_num=4;
        overlap_pixel=5;  % a half
       
        % padding A 
        A_pad=zeros(row+2*overlap_pixel, col+2*overlap_pixel, channel);
        A_pad(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=A;
        
        % do patch-wise completion
        for i=2:scale_num
            all_rr = 0;
            patch_size=[row/(2^(i-1))+2*overlap_pixel, col/(2^(i-1))+2*overlap_pixel, channel];
            conv_mat1=1/2*ones(patch_size);
            conv_mat1(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=ones(patch_size(1)-2*overlap_pixel,patch_size(2)-2*overlap_pixel,patch_size(3));
            conv_mat2=1/2*ones(patch_size);
            conv_mat2(overlap_pixel+1:end-overlap_pixel,overlap_pixel+1:end-overlap_pixel,:)=zeros(patch_size(1)-2*overlap_pixel,patch_size(2)-2*overlap_pixel,patch_size(3));
            
            num_patch = 0;
            for row_wise=1:2^(i-1)
                for col_wise=1:2^(i-1)
                    num_patch = num_patch + 1;
                    A_temp = A_pad((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:);
                    A_rr = relative_rank(A_temp);
                    all_rr = all_rr + A_rr;
                    % compare the patchs
                    figure(5);
                    imshow(A_pad((row_wise-1)*(patch_size(1)-2*overlap_pixel)+1:(row_wise-1)*(patch_size(1)-2*overlap_pixel)+patch_size(1),(col_wise-1)*(patch_size(2)-2*overlap_pixel)+1:(col_wise-1)*(patch_size(2)-2*overlap_pixel)+patch_size(2),:))
                    title('Original Patch');          
                end
            end
            patch_size = size(A_temp);
            fprintf('Number of patches in this stage is: %d. \n', num_patch)
            fprintf('Patch size is: %d \n', patch_size(1))
            avg_rr = all_rr / num_patch / patch_size(1);
            avg_rr_record{ff, i} = avg_rr;
        end
        fprintf('per image per miss rate costs time:   ')
        toc
    end
end