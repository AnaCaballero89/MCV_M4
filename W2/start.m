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
% 11 castle halloween pumpkins

indexPhoto = 10;

% version
% 1 initial
% 2 mixed gradients mean
% 3 mixed gradients laplacian (not work)

version = 2;

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
            destination = 'lena';
        case 4
            origin = 'nilai';
            destination = 'lena';
        case 5
            origin = 'girl';
            destination = 'ana';
        case 6
            origin = 'girl';
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
    
    
elseif indexPhoto == 9
    dst = double(imread('hojas.png'));
    src = double(imrotate(imread('mariquita.png'),90));
    % background
    aux1 = zeros(256,256);
    aux1(150:218, 100:245) = 1;
    mask_dst1 = logical(aux1);
    % ladybug
    aux1 = zeros(256,256);
    aux1(110:178, 85:230) = 1;
    mask_src1 = logical(aux1);
    mask_src1;
    
elseif indexPhoto == 10
    dst = double(imread('background1.png'));
    src = double(imread('toinsert1.png'));
    
    aux1 = zeros(256,256);
    aux1(188:255, 2:96) = 1;
    mask_dst1 = logical(aux1);
    % duck
    aux1 = zeros(256,256);
    aux1(45:112, 106:200) = 1;
    mask_src1 = logical(aux1);
    mask_src1;
else% indexPhoto == 11
    dst = double(imread('background2.png'));
    src = double(imread('toinsert2.png'));
    
    aux1 = zeros(256,256);
    aux1(15:175, 140:185) = 1;
    mask_dst1 = logical(aux1);
    % castle
    aux1 = zeros(256,256);
    aux1(15:175, 25:70) = 1;
    mask_src1 = logical(aux1);
    mask_src1;
end

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

%Preallocate
dst1 = src;
dst2 = src;

switch (version)
    case 1
        for nC = 1: nChannels
            
            drivingGrad_i = G5_DiFwd(src(:,:,nC),param.hi) - G5_DiBwd(src(:,:,nC),param.hi); % driving on i axis
            drivingGrad_j = G5_DjFwd(src(:,:,nC),param.hj) - G5_DjBwd(src(:,:,nC),param.hj); % driving on j axis
            
            driving_on_src = drivingGrad_i + drivingGrad_j; % union of the two gradients
            
            driving_on_dst = zeros(size(src(:,:,1)));
            driving_on_dst(mask_dst1) = driving_on_src(mask_src1);
            if indexPhoto < 9
                driving_on_dst(mask_dst2) = driving_on_src(mask_src2);
                mask_dst = mask_dst1 + mask_dst2;
            else
                mask_dst = mask_dst1;
            end
            
            param.driving = driving_on_dst;
            
            dst1(:,:,nC) = G5_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
            dst2(:,:,nC) = G5_Poisson_Equation_GaussSeidel(dst1(:,:,nC), mask_dst,  param);
        end
        
    case 2
        for nC = 1: nChannels
            % driving of the source image
            drivingGrad_i = G5_DiFwd(src(:,:,nC),param.hi) - G5_DiBwd(src(:,:,nC),param.hi); % driving on i axis
            drivingGrad_j = G5_DjFwd(src(:,:,nC),param.hj) - G5_DjBwd(src(:,:,nC),param.hj); % driving on j axis
            
            driving_on_src = drivingGrad_i + drivingGrad_j; % union of the two gradients
            
            driving_on_dst = zeros(size(src(:,:,1)));
            driving_on_dst(mask_dst1) = driving_on_src(mask_src1);
            if indexPhoto < 9
                driving_on_dst(mask_dst2) = driving_on_src(mask_src2);
            end
            
            %driving of the destination image
            drivingGrad_i = G5_DiFwd(dst(:,:,nC),param.hi) - G5_DiBwd(dst(:,:,nC),param.hi); % driving on i axis
            drivingGrad_j = G5_DjFwd(dst(:,:,nC),param.hj) - G5_DjBwd(dst(:,:,nC),param.hj); % driving on j axis
            
            driving_dst1 = drivingGrad_i + drivingGrad_j; % union of the two gradients
            
            driving_dst = zeros(size(dst(:,:,1)));
            driving_dst(mask_dst1) = driving_dst1(mask_dst1);
            if indexPhoto < 9
                driving_on_dst(mask_dst2) = driving_on_src(mask_src2);
                mask_dst = mask_dst1 + mask_dst2;
            else
                mask_dst = mask_dst1;
            end
            
            param.driving = (driving_on_dst + driving_dst)./2; % mean of the driving of source and destination, mixing gradients
            
            
            dst1(:,:,nC) = G5_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
            dst2(:,:,nC) = G5_Poisson_Equation_GaussSeidel(dst1(:,:,nC), mask_dst,  param);
        end
        
    case 3
        laplacian_dst = zeros(ni,nj);
        laplacian_on_dst = zeros(ni,nj);
        driving = zeros(ni,nj);
        mask_driving = zeros(ni,nj);
        for nC = 1: nChannels
            % driving of the source image
            drivingGrad_i = G5_DiFwd(src(:,:,nC),param.hi) - G5_DiBwd(src(:,:,nC),param.hi); % driving on i axis
            drivingGrad_j = G5_DjFwd(src(:,:,nC),param.hj) - G5_DjBwd(src(:,:,nC),param.hj); % driving on j axis
            
            driving_on_src = drivingGrad_i + drivingGrad_j; % union of the two gradients
            
            driving_on_dst = zeros(size(src(:,:,1)));
            driving_on_dst(mask_dst1) = driving_on_src(mask_src1);
            if indexPhoto < 9
                driving_on_dst(mask_dst2) = driving_on_src(mask_src2);
            end
            
            %driving of the destination image
            drivingGrad_i = G5_DiFwd(dst(:,:,nC),param.hi) - G5_DiBwd(dst(:,:,nC),param.hi); % driving on i axis
            drivingGrad_j = G5_DjFwd(dst(:,:,nC),param.hj) - G5_DjBwd(dst(:,:,nC),param.hj); % driving on j axis
            
            driving_dst1 = drivingGrad_i + drivingGrad_j; % union of the two gradients
            
            driving_dst = zeros(size(dst(:,:,1)));
            driving_dst(mask_dst1) = driving_dst1(mask_dst1);
            if indexPhoto < 9
                driving_dst(mask_dst2) = driving_dst1(mask_dst2);
                mask_dst = mask_dst1 + mask_dst2;
            else
                mask_dst = mask_dst1;
            end
            
            % laplacian of the gradients
            for i = 2:ni-1
                for j = 2:nj-1
                    if mask_dst1(i,j)
                        laplacian_on_dst(i,j) = driving_on_dst(i-1,j) + driving_on_dst(i+1,j) + driving_on_dst(i,j-1) + driving_on_dst(i,j+1) - 4*driving_on_dst(i,j);
                        laplacian_dst(i,j) = driving_dst(i-1,j) + driving_dst(i+1,j) + driving_dst(i,j-1) + driving_dst(i,j+1) - 4*driving_dst(i,j);
                        if (abs(laplacian_on_dst(i,j)) < abs(laplacian_dst(i,j)))
                            driving(i,j) = laplacian_dst(i,j);
                            mask_driving(i,j) = 1;
                        end
                    end
                end
            end
            
            param.driving = (laplacian_on_dst + laplacian_dst)./2;
            
            dst1(:,:,nC) = G5_Poisson_Equation_Axb(dst(:,:,nC), mask_driving,  param);
            dst2(:,:,nC) = G5_Poisson_Equation_GaussSeidel(dst1(:,:,nC), mask_driving,  param);
        end
end

figure, imshow(dst1/256);
figure, imshow(dst2/256);


