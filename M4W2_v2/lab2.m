%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 2: Image mosaics

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

imargb = imread('Data/llanes/llanes_a.jpg');
imbrgb = imread('Data/llanes/llanes_b.jpg');
imcrgb = imread('Data/llanes/llanes_c.jpg');

% imargb = imread('Data/castle_int/0016_s.png');
% imbrgb = imread('Data/castle_int/0015_s.png');
% imcrgb = imread('Data/castle_int/0014_s.png');

% imargb = imread('Data/aerial/site13/frame00000.png');
% imbrgb = imread('Data/aerial/site13/frame00002.png');
% imcrgb = imread('Data/aerial/site13/frame00003.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% imargb = double(imread('Data/aerial/site22/frame_00001.tif'));
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'));
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'));
% ima = imargb;
% imb = imbrgb;
% imc = imcrgb;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

figure;
imshow(imargb);%image(imargb)
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function --> done

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
    matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
H = eye(3);                             % --> done homography = 1 0 0; 0 1 0; 0 0 1
iwb = apply_H_v2(imbrgb, H, corners);   % ToDo: complete the call to the function 

H = Hab;                                % --> done Hab = ransac/homography2d(xab_a, xab_b)
iwa = apply_H_v2(imargb, H, corners);   % ToDo: complete the call to the function 

H = inv(Hbc);                           %  --> done Hbc = ransac/homography2d(xbc_b, xbc_c)
iwc = apply_H_v2(imcrgb, H, corners);   % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

% ToDo: compute the mosaic with castle_int images


% ToDo: compute the mosaic with aerial images set 13
% ToDo: compute the mosaic with aerial images set 22
% ToDo: comment the results in every of the four cases: say why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm
close all;
% Homography ab

x = points_a(1:2, matches_ab(1,inliers_ab));  % ToDo: set the non-homogeneous point coordinates of the --> done
xp = points_b(1:2, matches_ab(2,inliers_ab));  % point correspondences we will refine with the geometric method --> done
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

 Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 

err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);
 
Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm --> done 
% --> new part upt to: figure imshow(imargb)

xhat = P(10:end);
xhat = reshape(xhat, [2,size(xhat,1)/2]); 
% --> euclidean to homogeneous coordinates
xhat = [xhat ; ones(1,size(xhat,2))]; 
xhatp = Hab_r*xhat;

% --> euclidean coordinates from homogeneous
xh = xhat;
aux = size(xh, 1) - 1;
xhat = xh(1:aux,:) ./ repmat(xh(end,:), aux, 1);

% --> euclidean coordinates from homogeneous
xh = xhatp;
aux = size(xh, 1) - 1;
xhatp = xh(1:aux,:) ./ repmat(xh(end,:), aux, 1);



figure;
imshow(imargb);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
%% Homography bc
% ToDo: refine the homography bc with the Gold Standard algorithm
x = points_b(1:2, matches_bc(1,inliers_bc));  % ToDo: set the non-homogeneous point coordinates of the --> done
xp = points_c(1:2, matches_bc(2,inliers_bc));  % point correspondences we will refine with the geometric method --> done
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hbc(:) ; x(:) ];   

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 

err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%% Build mosaic

corners = [-400 1200 -100 650];
H = eye(3);                             
iwb = apply_H_v2(imbrgb, H, corners);   
H = Hab_r;                               
iwa = apply_H_v2(imargb, H, corners);   
H = inv(Hbc_r);                           
iwc = apply_H_v2(imcrgb, H, corners);  
figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
% %%              DLT algorithm (folder "logos").
% %%              Interpret and comment the results.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 7. OPTIONAL: Replace the logo of the UPF by the master logo
% %%              in one of the previous images using the DLT algorithm.



