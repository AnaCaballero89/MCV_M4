clearvars;

% indexPhoto
% 1 girl in Lena's face
% 2 Lena in girl's face
% 3 Ana in Lena's Face
% 4 Nilai in Lena's Face
% 5 Girl in Ana's Face
% 6 Girl in Nilai's Face
% 7 Ana in Nilai's Face
% 8 Nilai in Ana's Face
% 9 ladybug on leaves
% 10 duck and cascade

indexPhoto = 1;
if indexPhoto < 9
switch(indexPhoto) %choose name files for each index
    case 1
        origin = 'girl';
        destination = 'lena';
    case 2
        origin = 'lena';
        destination = 'girl';
    case 3
        origin = 'ana';
        destination = 'girl';
    case 4
        origin = 'nilai';
        destination = 'girl';
    case 5
        origin = 'lena';
        destination = 'ana';
    case 6
        origin = 'lena';
        destination = 'nilai';
    case 7
        origin = 'ana';
        destination = 'nilai';
    case 8
        origin = 'nilai';
        destination = 'ana';       
end

    src = double(imread(strcat(origin,'.png')));
    dst = double(imread(strcat(destination,'.png')));

    %masks to exchange: Eyes
    mask_src1 = logical(imread(strcat('mask_',origin,'_eyes.png')));
    mask_dst1 = logical(imread(strcat('mask_',destination,'_eyes.png')));

    %masks to exchange: Mouth
    mask_src2 = logical(imread(strcat('mask_',origin,'_mouth.png')));
    mask_dst2 = logical(imread(strcat('mask_',destination,'_mouth.png')));
   

elseif indexPhoto ==10
    dst = double(imread('hojas.png'));
    src = double(imread('mariquita.png'));
    aux1 = zeros(256,256);
    aux1(80:230, 80:150) = 1;
    mask_dst1 = logical(aux1);
    
    aux1 = zeros(256,256);
    aux1(80:230, 80:150) = 1;
    mask_src1 = logical(aux1);
    mask_src1;
    
else % indexPhoto==11
    dst = double(imread('background1.png'));
    src = double(imread('toinsert1.png'));
    
    aux1 = zeros(712,534);
    aux1(532:712, 1:161) = 1;
    mask_dst1 = logical(aux1);
    
    aux1 = zeros(712,534);
    aux1(120:300, 250:410) = 1;
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


