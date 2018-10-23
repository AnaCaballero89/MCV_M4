%Example script: You should replace the beginning of each function ('sol')
%with the name of your group. i.e. if your gropu name is 'G8' you should
%call :
% G8_DualTV_Inpainting_GD(I, mask, paramInp, paramROF)
close all;
clearvars;
%There are several black and white images to test:
%  image1_toRestore.jpg
%  image2_toRestore.jpg
%  image3_toRestore.jpg
%  image4_toRestore.jpg
%  image5_toRestore.jpg


% for to process images from 1 to 4
for imgNum = 1:4
    % ================================
    % ======   IMAGE =================
    % ================================

    %name= 'image5';
    %name= 'image1';
    name = strcat('image', int2str(imgNum));
    I = double(imread([ name '_toRestore.jpg']));
    %I=I(1:10,1:10);

    %Number of pixels for each dimension, and number of channles
    [ni, nj, nC] = size(I);

    if nC==3
        I = mean(I,3); %Convert to b/w. If you load a color image you should comment this line
    end

    %Normalize values into [0,1]
    I=I-min(I(:));
    I=I/max(I(:));

    % ================================
    % ======   MASK ==================
    % ================================

    %Load the mask
    mask_img = double(imread([name '_mask.jpg']));
    %mask_img =mask_img(1:10,1:10);
    [ni, nj, nC] = size(mask_img);
    if nC==3
        mask_img = mask_img(:,:,1); %Convert to b/w. If you load a color image you should comment this line
    end
    %We want to inpaint those areas in which mask == 1
    mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                          %mask(i,j) == 0 means we have information in that
                          %pixel


    % ================================
    % ====== GRADIENT DESCENT ========
    % ================================																	

    %%%Parameters for gradient descent (you do not need for week1)
    param.dt = 5*10^-7;
    param.iterMax = 10^4;
    param.tol = 10^-5;

    %%Parameters 
    param.hi = 1 / (ni-1);
    param.hj = 1 / (nj-1);

    % ================================
    % ====== LAPLACIAN + SHOW ========
    % ================================

    figure;%(1)
    subplot(1,2,1); imshow(I);
    title('Before')

    Iinp=Ana_Laplace_Equation_Axb(I, mask, param);
    %figure(2)
    subplot(1,2,2); imshow(Iinp);
    title('After');
    
    clearvars;
end
% ================================
% ================================
% Challenge image. (We have lost 99% of information)
 

%we make  plot to process images from 5 to 6
for imgNum = 5:6
    % ================================
    % ======   IMAGE =================
    % ================================
    if imgNum ==5
        img = 'image5_toRestore.jpg';
        msk = 'image5_mask.jpg';
    else
        img = 'image6_toRestore.tif';
        msk = 'image6_mask.tif';
    end
    
    I=double(imread(img)); 

    %Normalize values into [0,1]
    I=I/256;

    %Number of pixels for each dimension, and number of channels
    [ni, nj, nC] = size(I);

    % ================================
    % ======   MASK ==================
    % ================================

    mask_img=double(imread(msk));
    mask = mask_img >128; %mask(i,j) == 1 means we have lost information in that pixel
                          %mask(i,j) == 0 means we have information in that
                          %pixel

    % ================================
    % ====== GRADIENT DESCENT ========
    % ================================
    param.hi = 1 / (ni-1);
    param.hj = 1 / (nj-1);

    % ================================
    % ====== LAPLACIAN + SHOW ========
    % ================================
    figure;%(1)
    subplot(1,2,1);imshow(I);
    title('Before')

    Iinp(:,:,1)=Ana_Laplace_Equation_Axb(I(:,:,1), mask(:,:,1), param);
    Iinp(:,:,2)=Ana_Laplace_Equation_Axb(I(:,:,2), mask(:,:,2), param);
    Iinp(:,:,3)=Ana_Laplace_Equation_Axb(I(:,:,3), mask(:,:,3), param);

    %figure(2)
    subplot(1,2,2); imshow(Iinp);
    title('After');
    
    clearvars;
end
% ================================
% ================================
%% Goal Image
 clearvars;

 
% ================================
% ======   IMAGE =================
% ================================
%Read the image
I = double(imread('image_to_Restore.png'));

%Number of pixels for each dimension, and number of channels
[ni, nj, nC] = size(I);

%Normalize values into [0,1]
I = I - min(I(:));
I = I / max(I(:));

%We want to inpaint those areas in which mask == 1 (red part of the image)
I_ch1 = I(:,:,1);
I_ch2 = I(:,:,2);
I_ch3 = I(:,:,3);

% ================================
% ======   MASK ==================
% ================================
%TO COMPLETE 1
%mask_img(i,j) == 1 means we have lost information in that pixel
%mask(i,j) == 0 means we have information in that pixel
mask = zeros(size(I_ch1)); 
mask(I_ch1 == 1 & I_ch2 ==0 & I_ch3==0)=1; 


% ================================
% ====== GRADIENT DESCENT ========
% ================================
%%%Parameters for gradient descent (you do not need for week1)
%param.dt = 5*10^-7;
%param.iterMax = 10^4;
%param.tol = 10^-5;

%parameters
param.hi = 1 / (ni-1);
param.hj = 1 / (nj-1);

%% ================================
% ====== LAPLACIAN + SHOW ========
% ================================
figure;%(3)
subplot(1,2,1);
imshow(I);
title('Before')


Iinp(:,:,1)=Ana_Laplace_Equation_Axb(I_ch1, mask, param);
Iinp(:,:,2)=Ana_Laplace_Equation_Axb(I_ch2, mask, param);
Iinp(:,:,3)=Ana_Laplace_Equation_Axb(I_ch3, mask, param);

%figure(2)
subplot(1,2,2);
imshow(Iinp)
title('After');

