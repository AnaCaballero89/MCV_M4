clearvars;

% indexPhoto
% 1 example photo
% 2 example reverse
% 3 Ana in Lena's Face
% 4 Nilai in Lena's Face
% 5 Girl in Ana's Face
% 7 Girl in Nilai's Face
% 8 Ana in Nilai'Face
% 9 Nilai in Ana's Face
% 10 Duck in Waterfall (img1)
% 11 Squirrel in Garden (img2)

indexPhoto =1;
if indexPhoto ==1
    dst = double(imread('lena.png'));
    src = double(imread('girl.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    mask_src1=logical(imread('mask_src_eyes.png'));
    mask_dst1=logical(imread('mask_dst_eyes.png'));

    %masks to exchange: Mouth
    mask_src2=logical(imread('mask_src_mouth.png'));
    mask_dst2=logical(imread('mask_dst_mouth.png'));
elseif indexPhoto ==2
    dst = double(imread('girl.png'));
    src = double(imread('lena.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    mask_src1=logical(imread('mask_dst_eyes.png'));
    mask_dst1=logical(imread('mask_src_eyes.png'));

    %masks to exchange: Mouth
    mask_src2=logical(imread('mask_dst_mouth.png'));
    mask_dst2=logical(imread('mask_src_mouth.png'));
else % hay un problema con la mascara (parece)
    dst = double(imread('lena.png'));
    src = double(imread('ana.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    %mask_src1=logical(imread('mask_src_eyes.png'));
    %mask_dst1=logical(imread('mask_ana_eyes.png'));
 
    %masks to exchange: Mouth
    %mask_src2=logical(imread('mask_src_mouth.png'));
    %mask_dst2=logical(imread('mask_ana_mouth.png'));
    

end

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

%Preallocate
dst1= src;

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = G5_DiBwd(G5_DiFwd(src(:,:,nC),param.hi)); %??
    drivingGrad_j = G5_DjBwd(G5_DjFwd(src(:,:,nC),param.hj)); %??

    driving_on_src = drivingGrad_i + drivingGrad_j; %??
    
    driving_on_dst = zeros(size(src(:,:,1)));   
    driving_on_dst(mask_dst1(:)) = driving_on_src(mask_src1(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G5_Poisson_Equation_Axb(dst(:,:,nC), mask_dst1,  param);
end


for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = G5_DiBwd(G5_DiFwd(src(:,:,nC),param.hi)); %??
    drivingGrad_j = G5_DjBwd(G5_DjFwd(src(:,:,nC),param.hj)); %??

    driving_on_src = drivingGrad_i + drivingGrad_j; %??
    
    driving_on_dst = zeros(size(src(:,:,1)));  
    driving_on_dst(mask_dst2(:)) = driving_on_src(mask_src2(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G5_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst2,  param);
end

imshow(dst1/256)


