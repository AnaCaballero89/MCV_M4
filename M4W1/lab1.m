%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
% --> H structrue explained in slide 7 (Gloria's); H=[A t; v 1]
% --> H (3x3 Matrix) is a given homography 
% --> Rotation Matrix [ cos(theta) - sin(theta); cos(theta) sin(theta)]; 
% --> where theta is the orientation angle; Gloria's cass notes
s = 0.5;
theta = 30;
A = [(s*cosd(theta))  (-s*sind(theta)) ; (s*sind(theta))  (s*cosd(theta))];
t1 = 5;
t2 = 10;
v1 = 0;
v2 = 0;
% --> H = [A t ; 0 1]
% --> A is a non-singular 2 � 2 matrix
% --> T (t1 and t2) is a translation vector
% --> 0 (v1 and v2) is a 0 vector
H =  [A(1,1), A(1,2), t1 ; A(2,1), A(2,2), t2; v1, v2, 1];


I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
A = [0 1; 1 1];
t1 = 4;
t2 = 6;
% --> H = [A t ; 0 1]
% --> A is a non-singular 2 � 2 matrix
% --> T (t1 and t2) is a translation vector
% --> 0 (v1 and v2) is a 0 vector
H = [A(1,1), A(1,2), t1 ; A(2,1), A(2,2), t2; v1, v2, 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
% --> translation
T = [1 0 t1;0 1 t2; 0 0 1];

% --> rotation
% --> slide 6 lecture2a.pdf
[U,D,V] = svd(A);
Rtheta = U';
Rphi = V';
Rtheta = [Rtheta(1,1), Rtheta(1,2), 0 ; Rtheta(2,1), Rtheta(2,2), 0; 0,0,1];
Rphi = [Rphi(1,1), Rphi(1,2), 0 ; Rphi(2,1), Rphi(2,2), 0; 0,0,1];

% --> scale
D = [D(1,1), D(1,2), 0 ; D(2,1), D(2,2), 0 ; 0,0,1];

% H--> H decomposition
H_decomposition = T*(Rtheta*Rphi*D);

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
diff = round(H_decomposition)-H;

if (sum(diff(:))) == 0
    disp('matrices are equal')
else
    disp('matrices are not equal')
end
% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before

I2 = apply_H(I, H_decomposition);
figure; imshow(I); figure; imshow(uint8(I2));

diff = I2-I2_decomposition;
if diff == 0
    disp('images are equal')
else
    disp('images are not equal')
end

%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)


