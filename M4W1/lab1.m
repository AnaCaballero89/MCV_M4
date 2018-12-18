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
A = [(s*cosd(theta))  (-s*sind(theta)) ; 
    (s*sind(theta))  (s*cosd(theta))];
t1 = 5;
t2 = 10;
v1 = 0;
v2 = 0;
% --> H = [A t ; 0 1]
% --> A is a non-singular 2 � 2 matrix
% --> T (t1 and t2) is a translation vector
% --> 0 (v1 and v2) is a 0 vector
H =  [A(1,1),   A(1,2),     t1 ; 
      A(2,1),   A(2,2),     t2; 
      v1,       v2,          1];


I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
A = [0 1; 
    1 1];
t1 = 4;
t2 = 6;
% --> H = [A t ; 0 1]
% --> A is a non-singular 2 � 2 matrix
% --> T (t1 and t2) is a translation vector
% --> 0 (v1 and v2) is a 0 vector
H = [A(1,1),    A(1,2),     t1 ; 
    A(2,1),     A(2,2),     t2; 
    v1,         v2,          1];


I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
% --> translation
T = [1 0 t1;
    0 1 t2; 
    0 0 1];

% --> rotation
% --> slide 6 lecture2a.pdf
[U,D,V] = svd(A);
Rtheta = U*V';
Rphi = V';

Rtheta = [Rtheta(1,1),  Rtheta(1,2),    0 ; 
        Rtheta(2,1),    Rtheta(2,2),    0; 
        0,              0,              1];

Rphi = [Rphi(1,1),  Rphi(1,2), 0 ; 
        Rphi(2,1),   Rphi(2,2), 0; 
        0,          0,          1];

Rphit = Rphi';

% --> scale
S = [D(1,1),     D(1,2), 0 ; 
    D(2,1),      D(2,2), 0 ; 
    0,           0,      1];

% H--> H decomposition
H_decomposition = T*(Rtheta*Rphit*S*Rphi);

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

I2_decomposition = apply_H(I, H_decomposition);
figure; imshow(I); figure; imshow(uint8(I2_decomposition));

diff = I2-I2_decomposition;

if ((sum(diff(:))) == 0 ) || (isnan((sum(diff(:)))) ) 
    disp('images are equal')
else
    disp('images are not equal')
end

%% 1.3 Projective transformations (homographies)


% ToDo: generate a matrix H which produces a projective transformation,that
% can be descomposed in a product of 3 matrix, H= Hs*Ha*Hp, where Hs is a
% similarity transformation, Ha affinity transformation and Hp the
% projective one.

% First, the Hs similarity transformation, 
% Remember: Hs=[A t; 0 1], where A is the Rotation Matrix [ cos(theta) -
% sin(theta); cos(theta) sin(theta)]; and t the translation vector(t1,t2)
theta = 30;
s=0.5;
A = [(s*cosd(theta))  (-s*sind(theta)) ; (s*sind(theta))  (s*cosd(theta))];
t1 = 5;
t2 = 10;
v1 = 0;
v2 = 0;
Hs=  [A(1,1), A(1,2), t1 ; A(2,1), A(2,2), t2; v1, v2, 1];
% Then, the Ha affinity transformation 
% Remember: Ha=[K t ; 0 1], where K is a non singular 2x2 matrix and T (t1 
% and t2) is a translation vector 
K= [0 1; 
    1 1];
t1 = 0;
t2 = 0;
Ha= [K(1,1),    K(1,2),     t1 ; 
    K(2,1),     K(2,2),     t2; 
     0,         0,          1];
% Finally, we have Hp projective transformation. 
% Hp=[I 0; v vi], where v is not null 
Iden=[1 0;
   0 1];
v=[0 0];
vi=1;
Hp=[Iden [0;0]; v vi];

% Homography as H=Hs*Ha*Hp
H=Hs*Ha*Hp;

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

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different
% pairs of points

% Longer way

% coeff = polyfit([p1(1), p2(1)], [p1(2), p2(2)], 1); % y=mx+n 
% l1 = [-coeff(1) 1 -coeff(2)]; % (in images y increases downwards,y=-mx-n )
% coeff = polyfit([p3(1), p4(1)], [p3(2), p4(2)], 1); % y=mx+n 
% l2 = [-coeff(1) 1 -coeff(2)]; % (in images y increases downwards,y=-mx-n )
% coeff = polyfit([p5(1), p6(1)], [p5(2), p6(2)], 1); % y=mx+n 
% l3 = [-coeff(1) 1 -coeff(2)]; % (in images y increases downwards,y=-mx-n )
% coeff = polyfit([p7(1), p8(1)], [p7(2), p8(2)], 1); % y=mx+n 
% l4 = [-coeff(1) 1 -coeff(2)]; % (in images y increases downwards,y=-mx-n )

% Easier way
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
scatter(p1(1),p1(2),'r'); scatter(p2(1),p2(2),'r'); 
scatter(p3(1),p3(2),'r'); scatter(p4(1),p4(2),'r'); 
scatter(p5(1),p5(2),'r'); scatter(p6(1),p6(2),'r'); 
scatter(p7(1),p7(2),'r'); scatter(p8(1),p8(2),'r'); 
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

% Compute Vanishing points as cross product of "parallel" lines
van1 = cross(l1,l2);
van1 = van1/van1(3);
van2 = cross(l3, l4);
van2 = van2/van2(3);

% Compute vanishing line as cross product of vanishing points
l_va = cross(van1, van2);

% Matrix H 
H = [1               0               0;
     0               1               0;
     l_va(1)/l_va(3) l_va(2)/l_va(3) 1];

I3 = apply_H(I, H); 
figure; imshow(uint8(I3));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
HP1 = H*p1; HP2 = H*p2; HP3 = H*p3; HP4 = H*p4;
HP5 = H*p5; HP6 = H*p6; HP7 = H*p7; HP8 = H*p8;

% We normalise the transformed points 
HP1 = HP1/HP1(3); HP2 = HP2/HP2(3); HP3 = HP3/HP3(3); HP4 = HP4/HP4(3);
HP5 = HP5/HP5(3); HP6 = HP6/HP6(3); HP7 = HP7/HP7(3); HP8 = HP8/HP8(3);

% Make the transformed lines
lr1 = cross(HP1,HP2);
lr2 = cross(HP3,HP4);
lr3 = cross(HP5,HP6);
lr4 = cross(HP7,HP8);

% show the transformed lines in the transformed image
figure;imshow(uint8(I3));hold on;
scatter(HP1(1),HP1(2),'r'); scatter(HP2(1),HP2(2),'r'); 
scatter(HP3(1),HP3(2),'r'); scatter(HP4(1),HP4(2),'r'); 
scatter(HP5(1),HP5(2),'r'); scatter(HP6(1),HP6(2),'r'); 
scatter(HP7(1),HP7(2),'r'); scatter(HP8(1),HP8(2),'r'); 
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% Slopes of lines before affine transformation
slope_l1 = -l1(1)/l1(2); slope_l2 = -l2(1)/l2(2);
slope_l3 = -l3(1)/l3(2); slope_l4 = -l4(1)/l4(2);

% Slopes of lines after affine transformation
slope_lr1 = -lr1(1)/lr1(2); slope_lr2 = -lr2(1)/lr2(2);
slope_lr3 = -lr3(1)/lr3(2); slope_lr4 = -lr4(1)/lr4(2);

% Angle between lines before and after transformation 
% (after is 0 as those lines are now parallel)
ang_l1_l2 = atand((slope_l1-slope_l2)/(1+slope_l1*slope_l2));
ang_lr1_lr2 = atand((slope_lr1-slope_lr2)/(1+slope_lr1*slope_lr2));

ang_l3_l4 = atand((slope_l3-slope_l4)/(1+slope_l3*slope_l4));
ang_lr3_lr4 = atand((slope_lr3-slope_lr4)/(1+slope_lr3*slope_lr4));

fprintf('Angle between l1 and l2 is , %f degrees\n', ang_l1_l2);
fprintf('Angle between affine rectified lr1 and lr2 is , %f degrees\n', ang_lr1_lr2);
fprintf('Angle between l3 and l4 is , %f degrees\n', ang_l3_l4);
fprintf('Angle between affine rectified lr3 and lr4 is , %f degrees\n', ang_lr3_lr4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification
%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rec�tification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.
% 
% --> Orthogonal pair of lines 1
l1 = lr1;
m1 = lr3;
% --> Orthogonal pair of lines 2
l2 = lr2;
m2 = lr4;

% --> slide 48 /lecture2_CB.pdf
A = [l1(1)*m1(1),   l1(1)*m1(2)+l1(2)*m1(1),    l1(2)*m1(2);
     l2(1)*m2(1),   l2(1)*m2(2)+l2(2)*m2(1),    l2(2)*m2(2)];
 % --> orthonormal basis for the null space of A 
s_vec = null(A);
% --> solve a system of equations to get S
S = [s_vec(1),  s_vec(2); 
    s_vec(2),   s_vec(3)];

K = chol(S); % --> upper triangular matrix K from the diagonal and upper triangle of matrix S
H = eye(3); % --> 3x3 identity matrix with ones on the main diagonal and zeros elsewhere
K = inv(K); % --> inverse matrix
H(1:2,1:2) = K;
H = H';

% --> transformed lines
l1trans = H'\l1;%inv(H')*l1;
l2trans = H'\l2;%inv(H')*l2;
m1trans = H'\m1;%inv(H')*m1;
m2trans = H'\m2;%inv(H')*m2;

% --> VISUALIZE LINES
% --> visualize original lines
figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t+l1(3))/l1(2), 'y');
plot(t, -(m1(1)*t+m1(3))/m1(2), 'y');
plot(t, -(l2(1)*t+l2(3))/l2(2), 'y');
plot(t, -(m2(1)*t+m2(3))/m2(2), 'y');

% --> visualize transformed lines
I4 = apply_H(I3, H); % --> I3 comes from previous part
figure; imshow(uint8(I4));
hold on;
t=1:0.1:1000;
plot(t, -(l1trans(1)*t+l1trans(3))/l1trans(2), 'y');
plot(t, -(m1trans(1)*t+m1trans(3))/m1trans(2), 'y');
plot(t, -(l2trans(1)*t+l2trans(3))/l2trans(2), 'y');
plot(t, -(m2trans(1)*t+m2trans(3))/m2trans(2), 'y');

% Compute angles Test --------------begin--------------------
% Slopes of lines before affine transformation
slope_l1 = -l1(1)/l1(2); slope_l2 = -l2(1)/l2(2);
slope_m1 = -m1(1)/m1(2); slope_m2 = -m2(1)/m2(2);

% Slopes of lines after affine transformation
slope_l1trans = -l1trans(1)/l1trans(2); slope_l2trans = -l2trans(1)/l2trans(2);
slope_m1trans = -m1trans(1)/m1trans(2); slope_m2trans = -m2trans(1)/m2trans(2);

% Angle between lines before and after transformation 
ang_l1_m1 = atand((slope_l1-slope_m1)/(1+slope_l1*slope_m1));
ang_l1trans_m1trans = atand((slope_l1trans-slope_m1trans)/(1+slope_l1trans*slope_m1trans));

ang_l2_m2 = atand((slope_l2-slope_m2)/(1+slope_l2*slope_m2));
ang_l2trans_m2trans = atand((slope_l2trans-slope_m2trans)/(1+slope_l2trans*slope_m2trans));

fprintf('Angle between l1 and m1 is , %f degrees\n', ang_l1_m1);
fprintf('Angle between affine+metric rectified l1 transf. and m1 transf. is , %f degrees\n', ang_l1trans_m1trans);
fprintf('Angle between l2 and m2 is , %f degrees\n', ang_l2_m2);
fprintf('Angle between affine+metric rectified l2 transf. and m2 transf. is , %f degrees\n', ang_l2trans_m2trans);
% Test --------------end--------------------


% % --> ORIGINAL LINES
% % --> Normalize
% l1 = l1/l1(3);
% m1 = m1/m1(3);
% degrees = atan2d(norm(cross(l1,m1)),dot(l1,m1));
% %degrees = angle(l1(1:2), m1(1:2));
% fprintf('angle between metric rectified l1 and m1 is %f degrees \n', degrees);
% 
% % --> Normalize
% l2 = l2/l2(3);
% m2 = m2/m2(3);
% degrees = atan2d(norm(cross(l2,m2)),dot(l2,m2));
% %degrees = angle(l2(1:2), m2(1:2));
% fprintf('angle between metric rectified l2 and m2 is %f degrees \n', degrees);
% 
% % --> TRANSOFRMED LINES 
% % --> Normalize
% l1trans = l1trans/l1trans(3);
% m1trans = m1trans/m1trans(3);
% degrees = atan2d(norm(cross(l1trans,m1trans)),dot(l1trans,m1trans));
% %degrees = angle(l1trans(1:2), m1trans(1:2));
% fprintf('angle between metric rectified l1trans and m1trans is %f degrees \n', degrees);
% 
% % --> Normalize
% l2trans = l2trans/l2trans(3);
% m2trans = m2trans/m2trans(3);
% degrees = atan2d(norm(cross(l2trans,m2trans)),dot(l2trans,m2trans));
% %degrees = angle(l2trans(1:2), m2trans(1:2));
% fprintf('angle between metric rectified l2trans and m2trans is %f degrees \n', degrees);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.
close all;
I=imread('Data/0001_s_crop.png');
A = load('Data/0001_s_info_lines.txt');

% indices of lines
i = 159;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 614;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 541;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 645;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% Easier way
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
scatter(p1(1),p1(2),'r'); scatter(p2(1),p2(2),'r'); 
scatter(p3(1),p3(2),'r'); scatter(p4(1),p4(2),'r'); 
scatter(p5(1),p5(2),'r'); scatter(p6(1),p6(2),'r'); 
scatter(p7(1),p7(2),'r'); scatter(p8(1),p8(2),'r'); 
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

% Compute Vanishing points as cross product of "parallel" lines
van1 = cross(l1,l2);
van1 = van1/van1(3);
van2 = cross(l3, l4);
van2 = van2/van2(3);

% Compute vanishing line as cross product of vanishing points
l_va = cross(van1, van2);

% Matrix H 
H = [1               0               0;
     0               1               0;
     l_va(1)/l_va(3) l_va(2)/l_va(3) 1];

I3 = apply_H(I, H); 
figure; imshow(uint8(I3));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
HP1 = H*p1; HP2 = H*p2; HP3 = H*p3; HP4 = H*p4;
HP5 = H*p5; HP6 = H*p6; HP7 = H*p7; HP8 = H*p8;

% We normalise the transformed points 
HP1 = HP1/HP1(3); HP2 = HP2/HP2(3); HP3 = HP3/HP3(3); HP4 = HP4/HP4(3);
HP5 = HP5/HP5(3); HP6 = HP6/HP6(3); HP7 = HP7/HP7(3); HP8 = HP8/HP8(3);

% Make the transformed lines
lr1 = cross(HP1,HP2);
lr2 = cross(HP3,HP4);
lr3 = cross(HP5,HP6);
lr4 = cross(HP7,HP8);

% show the transformed lines in the transformed image
figure;imshow(uint8(I3));hold on;
scatter(HP1(1),HP1(2),'r'); scatter(HP2(1),HP2(2),'r'); 
scatter(HP3(1),HP3(2),'r'); scatter(HP4(1),HP4(2),'r'); 
scatter(HP5(1),HP5(2),'r'); scatter(HP6(1),HP6(2),'r'); 
scatter(HP7(1),HP7(2),'r'); scatter(HP8(1),HP8(2),'r'); 
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% Slopes of lines before affine transformation
slope_l1 = -l1(1)/l1(2); slope_l2 = -l2(1)/l2(2);
slope_l3 = -l3(1)/l3(2); slope_l4 = -l4(1)/l4(2);

% Slopes of lines after affine transformation
slope_lr1 = -lr1(1)/lr1(2); slope_lr2 = -lr2(1)/lr2(2);
slope_lr3 = -lr3(1)/lr3(2); slope_lr4 = -lr4(1)/lr4(2);

% Angle between lines before and after transformation 
% (after is 0 as those lines are now parallel)
ang_l1_l2 = atand((slope_l1-slope_l2)/(1+slope_l1*slope_l2));
ang_lr1_lr2 = atand((slope_lr1-slope_lr2)/(1+slope_lr1*slope_lr2));

ang_l3_l4 = atand((slope_l3-slope_l4)/(1+slope_l3*slope_l4));
ang_lr3_lr4 = atand((slope_lr3-slope_lr4)/(1+slope_lr3*slope_lr4));

fprintf('Angle between l1 and l2 is , %f degrees\n', ang_l1_l2);
fprintf('Angle between affine rectified lr1 and lr2 is , %f degrees\n', ang_lr1_lr2);
fprintf('Angle between l3 and l4 is , %f degrees\n', ang_l3_l4);
fprintf('Angle between affine rectified lr3 and lr4 is , %f degrees\n', ang_lr3_lr4);

% Metric rectification
l1 = lr1;
m1 = lr3;
% --> Orthogonal pair of lines 2
l2 = lr2;
m2 = lr4;

% --> slide 48 /lecture2_CB.pdf
A = [l1(1)*m1(1),   l1(1)*m1(2)+l1(2)*m1(1),    l1(2)*m1(2);
     l2(1)*m2(1),   l2(1)*m2(2)+l2(2)*m2(1),    l2(2)*m2(2)];
 % --> orthonormal basis for the null space of A 
s_vec = null(A);
% --> solve a system of equations to get S
S = [s_vec(1),  s_vec(2); 
    s_vec(2),   s_vec(3)];

K = chol(S); % --> upper triangular matrix K from the diagonal and upper triangle of matrix S
H = eye(3); % --> 3x3 identity matrix with ones on the main diagonal and zeros elsewhere
K = inv(K); % --> inverse matrix
H(1:2,1:2) = K;
H = H';

% --> transformed lines
l1trans = H'\l1;%inv(H')*l1;
l2trans = H'\l2;%inv(H')*l2;
m1trans = H'\m1;%inv(H')*m1;
m2trans = H'\m2;%inv(H')*m2;


% --> VISUALIZE LINES
% --> visualize original lines
figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t+l1(3))/l1(2), 'y');
plot(t, -(m1(1)*t+m1(3))/m1(2), 'y');
plot(t, -(l2(1)*t+l2(3))/l2(2), 'y');
plot(t, -(m2(1)*t+m2(3))/m2(2), 'y');

% --> visualize transformed lines
I4 = apply_H(I3, H); % --> I3 comes from previous part
figure; imshow(uint8(I4));
hold on;
t=1:0.1:1000;
plot(t, -(l1trans(1)*t+l1trans(3))/l1trans(2), 'y');
plot(t, -(m1trans(1)*t+m1trans(3))/m1trans(2), 'y');
plot(t, -(l2trans(1)*t+l2trans(3))/l2trans(2), 'y');
plot(t, -(m2trans(1)*t+m2trans(3))/m2trans(2), 'y');

% Compute angles Test --------------begin--------------------
% Slopes of lines before affine transformation
slope_l1 = -l1(1)/l1(2); slope_l2 = -l2(1)/l2(2);
slope_m1 = -m1(1)/m1(2); slope_m2 = -m2(1)/m2(2);

% Slopes of lines after affine transformation
slope_l1trans = -l1trans(1)/l1trans(2); slope_l2trans = -l2trans(1)/l2trans(2);
slope_m1trans = -m1trans(1)/m1trans(2); slope_m2trans = -m2trans(1)/m2trans(2);

% Angle between lines before and after transformation 
ang_l1_m1 = atand((slope_l1-slope_m1)/(1+slope_l1*slope_m1));
ang_l1trans_m1trans = atand((slope_l1trans-slope_m1trans)/(1+slope_l1trans*slope_m1trans));

ang_l2_m2 = atand((slope_l2-slope_m2)/(1+slope_l2*slope_m2));
ang_l2trans_m2trans = atand((slope_l2trans-slope_m2trans)/(1+slope_l2trans*slope_m2trans));

fprintf('Angle between l1 and m1 is , %f degrees\n', ang_l1_m1);
fprintf('Angle between affine+metric rectified l1 transf. and m1 transf. is , %f degrees\n', ang_l1trans_m1trans);
fprintf('Angle between l2 and m2 is , %f degrees\n', ang_l2_m2);
fprintf('Angle between affine+metric rectified l2 transf. and m2 transf. is , %f degrees\n', ang_l2trans_m2trans);
% Test --------------end--------------------


% % --> ORIGINAL LINES
% % --> Normalize
% l1 = l1/l1(3);
% m1 = m1/m1(3);
% degrees = atan2d(norm(cross(l1,m1)),dot(l1,m1));
% %degrees = angle(l1(1:2), m1(1:2));
% fprintf('angle between metric rectified l1 and m1 is %f degrees \n', degrees);
% 
% % --> Normalize
% l2 = l2/l2(3);
% m2 = m2/m2(3);
% degrees = atan2d(norm(cross(l2,m2)),dot(l2,m2));
% %degrees = angle(l2(1:2), m2(1:2));
% fprintf('angle between metric rectified l2 and m2 is %f degrees \n', degrees);
% 
% % --> TRANSOFRMED LINES 
% % --> Normalize
% l1trans = l1trans/l1trans(3);
% m1trans = m1trans/m1trans(3);
% degrees = atan2d(norm(cross(l1trans,m1trans)),dot(l1trans,m1trans));
% %degrees = angle(l1trans(1:2), m1trans(1:2));
% fprintf('angle between metric rectified l1trans and m1trans is %f degrees \n', degrees);
% 
% % --> Normalize
% l2trans = l2trans/l2trans(3);
% m2trans = m2trans/m2trans(3);
% degrees = atan2d(norm(cross(l2trans,m2trans)),dot(l2trans,m2trans));
% %degrees = angle(l2trans(1:2), m2trans(1:2));
% fprintf('angle between metric rectified l2trans and m2trans is %f degrees \n', degrees);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



