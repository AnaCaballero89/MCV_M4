
clear all;
 
%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x2 = points{i}(1:2, matches(2, :));
    H{i} = 0;
    x1_m = [x1; ones(1, size(x1, 2))];
    x2_m = [x2; ones(1, size(x2, 2))];
    th=3;
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1_m, x2_m, th, 1000);

    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
    %vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
% 

for i = 1:N
    h = H{i};
    A(2*i - 1,:) = [h(1,1) * h(1,2),  h(1,1) * h(2,2) + h(2,1) * h(1,2), ...
                    h(1,1) * h(3,2) + h(3,1) * h(1,2), h(2,1) * h(2,2), ...
                    h(2,1) * h(3,2) + h(3,1) * h(2,2), h(3,1) * h(3,2)];
    A(2*i,:) = [h(1,1) * h(1,1),  h(1,1) * h(2,1) + h(2,1) * h(1,1), ...
                h(1,1) * h(3,1) + h(3,1) * h(1,1) , ...
                h(2,1) * h(2,1), h(2,1) * h(3,1) + h(3,1) * h(2,1), ...
                h(3,1) * h(3,1)]   -  [h(1,2) * h(1,2), ...
                h(1,2) * h(2,2) + h(2,2) * h(1,2), ...
                h(1,2) * h(3,2) + h(3,2) * h(1,2), ...
                h(2,2) * h(2,2), ...
                h(2,2) * h(3,2) + h(3,2) * h(2,2), ...
                h(3,2) * h(3,2)];
end
[U,S,V]=svd(A);
%w symmetric matrix, the solution of the system of equations is the last
%column of V so we get this:
w=[V(1,end),V(2,end),V(3,end);
   V(2,end),V(4,end),V(5,end);
   V(3,end),V(5,end),V(6,end)];
 
%% Recover the camera calibration.

K = inv(chol(w)); 
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
    h = H{i};
    r1 = inv(K)*h(:,1)/sqrt(sum((inv(K)*h(:,1)).^2 ));
    r2 = inv(K)*h(:,2)/sqrt(sum((inv(K)*h(:,2)).^2 ));
    t{i} = inv(K)*h(:,3)/sqrt(sum((inv(K)*h(:,1)).^2));
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
for i = 1:N
    vgg_scatter_plot( R{i}*[p1+t{i} p2+t{i} p3+t{i} p4+t{i} p1+t{i}], 'r');
    %vgg_scatter_plot(R{i}.*[p1 p2 p3 p4 zeros(1,5)], 'r');
end
% 
%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    X_m=[X; ones(1, size(X,2))];
    x = euclid(P{i} * (X_m));
    vgg_scatter_plot(x, 'g');
end
