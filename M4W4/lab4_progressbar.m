%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 4: Reconstruction from two views (knowing internal camera parameters) 

close all;
clc;

addpath('../M4W2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Triangulation
basicwaitbar = waitbar(0,'Please wait...');
pause(.5)

disp('... Runing exercice 1')
% ToDo: create the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

%% Test the triangulate function
% Use this code to validate that the function triangulate works properly
waitbar(0,basicwaitbar,'Ex1: Loading your data');
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;

X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

waitbar(0.08,basicwaitbar,'Ex1: Triangulating');
N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

waitbar(0.74,basicwaitbar,'Ex1: Calculating Euclidian Error');
error = euclid(X_test) - euclid(X_trian)

close(basicwaitbar');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reconstruction from two views
basicwaitbar = waitbar(0,'Please wait...');
pause(.5)

disp('... Runing exercice 2')
waitbar(0,basicwaitbar,'Ex2: Loading your data');
%% Read images
Irgb{1} = imread('Data/0001_s.png');
Irgb{2} = imread('Data/0002_s.png');
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});

waitbar(0.48,basicwaitbar,'Ex2: Computing key points');
disp('computing key points...')
%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

waitbar(0.51,basicwaitbar,'Ex2: ploting matches')
disp('ploting matches...')
% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');

waitbar(05648,basicwaitbar,'Ex2: fit fundamental matrix')
disp('fit fundamental matrix...')
%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

disp('plot inliers...')
% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

disp('computing candidate camera matrices...')
%% Compute candidate camera matrices.
% Camera calibration matrix
K = [2362.12 0 1520.69; 0 2366.12 1006.81; 0 0 1];
scale = 0.3;
H = [scale 0 0; 0 scale 0; 0 0 1];
K = H * K;

% ToDo: Compute the Essential matrix from the Fundamental matrix
E = K'*F*K; 

% ToDo: write the camera projection matrix for the first camera
P1 =K*[eye(3),zeros(3,1)]; 

%
[U,S,V] = svd(E);
W = [0 -1 0; 1 0 0; 0 0 1];
u3=U(:,end);

waitbar(0.59,basicwaitbar,'Ex2: writting 4 possible matrices')
disp('writting 4 possible matrices...')
% ToDo: write the four possible matrices for the second camera
Pc2 = {};
Pc2{1} =K*[U*W*V',u3]; 
Pc2{2} = K*[U*W*V',-u3];
Pc2{3} = K*[U*W'*V',u3];
Pc2{4} = K*[U*W'*V',-u3];

% HINT: You may get improper rotations; in that case you need to change
%       their sign.
% Let R be a rotation matrix, you may check:
% if det(R) < 0
%     R = -R;
% end
waitbar(0.94,basicwaitbar,'Ex2: plot 1st camera and solutions for 2nd camera')
disp('ploting 1st camera and solutions for 2nd camera...')
% plot the first camera and the four possible solutions for the second
figure;
plot_camera(P1,w,h);
plot_camera(Pc2{1},w,h);
plot_camera(Pc2{2},w,h);
plot_camera(Pc2{3},w,h);
plot_camera(Pc2{4},w,h);

%% Reconstruct structure

disp('Reconstructing structure...')
% ToDo: Choose a second camera candidate by triangulating a match.
for i=1:4
    P2 = Pc2{i}; % Pick a matching point
    trian = triangulate(x1(:,1), x2(:,1), P1, P2, [w h]); % triangulate the matching point
    proj1 = P1*trian; % project to the camera
    proj2 = P2*trian; % project to the camera
    if (proj1(3) >= 0) && (proj2(3) >= 0)
        correct3D = i; % correct camera matrix 
    end
end
P2 = Pc2{correct3D};

% Triangulate all matches.
waitbar(0.97,basicwaitbar,'Ex2: triangulating all matches')
disp('Triangulating all matches...')
N = size(x1,2);
X = zeros(4,N);
for i = 1:N
    X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
end

%% Plot with colors
waitbar(0.99,basicwaitbar,'Ex2: Plot with color')
disp('Plot with colors...')
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
plot_camera(P1,w,h);
plot_camera(P2,w,h);
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end;
axis equal;


%% Compute reprojection error.

disp('Computing reprojection error...')
% ToDo: compute the reprojection errors
projx1 = euclid(P1*X);
projx2 = euclid(P2*X);


disp('ploting histogram...')
%       plot the histogram of reprojection errors, and
projErrors1 = sqrt(sum((x1-projx1).^2, 1));
projErrors2 = sqrt(sum((x2-projx2).^2, 1));

histfit([projErrors1 projErrors2]);
hold on
totalErrorProjected_1 = sum(projErrors1)
totalErrorProjected_2 = sum(projErrors2)

%       plot the mean reprojection error
total_errorProjected = totalErrorProjected_1+totalErrorProjected_2;
n_points = size(x1,2);

disp('calculating mean error...')
meanError = (total_errorProjected/(n_points*2))
line([meanError meanError], ylim, 'Color','r');

close(basicwaitbar');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Depth map computation with local methods (SSD)
basicwaitbar = waitbar(0,'Please wait...');
pause(.5)

disp('... Runing exercice 3')
% Data images: 'scene1.row3.col3.ppm','scene1.row3.col4.ppm'
% Disparity ground truth: 'truedisp.row3.col3.pgm'

% Write a function called 'stereo_computation' that computes the disparity
% between a pair of rectified images using a local method based on a matching cost 
% between two local windows.
% 
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
%
% In this part we ask to implement only the SSD cost
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.
%
% Note 1: Use grayscale images
% Note 2: For this first set of images use 0 as minimum disparity 
% and 16 as the the maximum one.
waitbar(0,basicwaitbar,'Ex3: Loading your data');
disp('setup variables...')
leftImg = rgb2gray(imread('Data/scene1.row3.col3.ppm'));
rightImg = rgb2gray(imread('Data/scene1.row3.col4.ppm'));
%gt = imread('Data/truedisp.row3.col3.pgm');
minDisp = 0;  % minimum disparity; Note 1
maxDisp = 16; % maximum disparity; Note 2
mc = 'SSD';   % matching cost

waitbar(0,basicwaitbar,'Ex3: stereo computation, SSD for window size = 3x3'); 
disp('stereo_computation, SSD for window size = 3x3...')
ws = 3;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.15,basicwaitbar,'Ex3: stereo computation, SSD for window size = 9x9'); 
disp('stereo_computation, SSD for window size = 9x9...')
ws = 9;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.34,basicwaitbar,'Ex3: stereo computation, SSD for window size = 20x20'); 
disp('stereo_computation, SSD for window size = 20x20...')
ws = 20;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.61,basicwaitbar,'Ex3: stereo computation, SSD for window size = 30x30'); 
disp('stereo_computation, SSD for window size = 30x30...')
ws = 30;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

close(basicwaitbar');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 4. Depth map computation with local methods (NCC)
basicwaitbar = waitbar(0,'Please wait...');
pause(.5)
disp('... Runing exercice 4') 
% Complete the previous function by adding the implementation of the NCC
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.

% leftImg = rgb2gray(imread('Data/scene1.row3.col3.ppm'));
% rightImg = rgb2gray(imread('Data/scene1.row3.col4.ppm'));
% gt = imread('Data/truedisp.row3.col3.pgm');
% min_disp = 0;
% max_disp = 16;
waitbar(0,basicwaitbar,'Ex4: Loading your data');
mc = 'NCC';

waitbar(0,basicwaitbar,'Ex4: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 3x3...')
ws = 3;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.15,basicwaitbar,'Ex4: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 9x9...')
ws = 9;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.34,basicwaitbar,'Ex4: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 20x20...')
ws = 20;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.61,basicwaitbar,'Ex4: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 30x30...')
ws = 30;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max 

close(basicwaitbar');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Depth map computation with local methods
basicwaitbar = waitbar(0,'Please wait...');
pause(.5)
disp('... Runing exercice 5') 

% Data images: '0001_rectified_s.png','0002_rectified_s.png'

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.
waitbar(0,basicwaitbar,'Ex5: Loading your data');
disp('setup variables...')

leftImg = rgb2gray(imread('Data/0001_rectified_s.png'));
rightImg = rgb2gray(imread('Data/0002_rectified_s.png'));
minDisp = 0;  % minimum disparity; Note 1
maxDisp = 16; % maximum disparity; Note 2

mc = 'SSD';   % matching cost

waitbar(0.05,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, SSD for window size = 3x3...')
ws = 3;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.12,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, SSD for window size = 9x9...')
ws = 9;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.25,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, SSD for window size = 20x20...')
ws = 20;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.37,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, SSD for window size = 30x30...')
ws = 30;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity


mc = 'NCC';
waitbar(0.50,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 3x3...')
ws = 3;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.63,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 9x9...')
ws = 9;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.75,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 20x20...')
ws = 20;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

waitbar(0.88,basicwaitbar,'Ex5: stereo computation, SSD for window size = 3x3');
disp('stereo_computation, NCC for window size = 30x30...')
ws = 30;      % window size 
disparity = stereo_computation(leftImg, rightImg, minDisp, maxDisp, ws, mc);
figure; imshow(uint8(disparity)*16); % 16 'cos of the max disparity

close(basicwaitbar');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 6. Bilateral weights
% 
% % Modify the 'stereo_computation' so that you can use bilateral weights (or
% % adaptive support weights) in the matching cost of two windows.
% % Reference paper: Yoon and Kweon, "Adaptive Support-Weight Approach for Correspondence Search", IEEE PAMI 2006
% %
% % Comment the results and compare them to the previous results (no weights).
% %
% % Note: Use grayscale images (the paper uses color images)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OPTIONAL:  7. Stereo computation with Belief Propagation
% 
% % Use the UGM library used in module 2 and implement a  
% % stereo computation method that minimizes a simple stereo energy with 
% % belief propagation. 
% % For example, use an L2 or L1 pixel-based data term (SSD or SAD) and 
% % the same regularization term you used in module 2. 
% % Or pick a stereo paper (based on belief propagation) from the literature 
% % and implement it. Pick a simple method or just simplify the method they propose.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OPTIONAL:  8. Depth computation with Plane Sweeping
% 
% % Implement the plane sweeping method explained in class.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OPTIONAL:  9. Depth map fusion 
% 
% % In this task you are asked to implement the depth map fusion method
% % presented in the following paper:
% % B. Curless and M. Levoy. A Volumetric Method for Building Complex
% % Models from Range Images. In Proc. SIGGRAPH, 1996.
% %
% % 1. Use the set of facade images 00xx_s.png to compute depth maps 
% % corresponding to different views (and optionally from different pairs of 
% % images for the same view).
% % 2. Then convert each depth map to a signed distance function defined in 
% % a disretized volume (using voxels).
% % 3. Average the different signed distance functions, the resulting 
% % signed distance is called D.
% % 4. Set as occupied voxels (those representing the surface) those 
% % where D is very close to zero. The rest of voxels will be considered as 
% % empty.
% %
% % For that you need to compute a depth map from a pair of views in general
% % position (non rectified). Thus, you may either use the plane sweep
% % algorithm (if you did it) or the local method for estimating depth
% % (mandatory task) together with the following rectification method which 
% % has an online demo available: 
% % http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OPTIONAL:  10. New view synthesis
% 
% % In this task you are asked to implement part of the new view synthesis method
% % presented in the following paper:
% % S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.
% 
% % You will use a pair of rectified stereo images (no need for prewarping
% % and postwarping stages) and their corresponding ground truth disparities
% % (folder "new_view").
% % Remember to take into account occlusions as explained in the lab session.
% % Once done you can apply the code to the another pair of rectified images 
% % provided in the material and use the estimated disparities with previous methods.
