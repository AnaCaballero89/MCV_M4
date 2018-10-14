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
% 10 ladybug on leaves

indexPhoto =11;
if indexPhoto ==1 % <- works!
    dst = double(imread('lena.png'));
    src = double(imread('girl.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    mask_src1=logical(imread('mask_src_eyes.png'));
    mask_dst1=logical(imread('mask_dst_eyes.png'));

    %masks to exchange: Mouth
    mask_src2=logical(imread('mask_src_mouth.png'));
    mask_dst2=logical(imread('mask_dst_mouth.png'));
elseif indexPhoto ==2% <- works!
    dst = double(imread('girl.png'));
    src = double(imread('lena.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    mask_src1=logical(imread('mask_dst_eyes.png'));
    mask_dst1=logical(imread('mask_src_eyes.png'));

    %masks to exchange: Mouth
    mask_src2=logical(imread('mask_dst_mouth.png'));
    mask_dst2=logical(imread('mask_src_mouth.png'));
elseif indexPhoto ==3% <- hay un problema con la mascara (parece),
    dst = double(imread('lena.png'));
    src = double(imread('ana.png')); % flipped girl, because of the eyes

    %masks to exchange: Eyes
    %mask_src1=logical(imread('mask_src_eyes.png'));
    mask_dst1=logical(imread('mask_dst_eyes.png'));

    %masks to exchange: Mouth
    %mask_src2=logical(imread('mask_src_mouth.png'));
    mask_dst2=logical(imread('mask_dst_mouth.png'));
    
    %%generate masks
    [a,b] =size(src);
    aux1 = zeros(a, b);
    aux1(256:351, 163:421) = 1;
    mask_src1 = logical(aux1);

    aux2 = zeros(a, b);
    aux2(417:479, 231:348) = 1;
    mask_src2 = logical(aux2);
elseif indexPhoto ==11% <- works!
    dst = double(imread('hojas.png'));
    src = double(imread('mariquita.png'));
    aux1 = zeros(256,256);
    aux1(80:230, 80:150) = 1;
    mask_dst1 = logical(aux1);
    
    aux1 = zeros(256,256);
    aux1(80:230, 80:150) = 1;
    mask_src1 = logical(aux1);
    mask_src1;
end

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

%Preallocate
dst1= src;
if indexPhoto <10 
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
else
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
end
imshow(dst1/256)


