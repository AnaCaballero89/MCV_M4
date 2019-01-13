function [ F ] = fundamental_matrix(x1_test, x2_test)
% Computation of the Fundamental Matrix with 8 points or more.
%   Detailed explanation goes here

numb=size(x1_test,2); 

%normalization of the points by the function provided 

[x1_norm, T1] = normalise2dpts(x1_test);
[x2_norm, T2] = normalise2dpts(x2_test);

%Compute linear combination of the matches such as W*f=0
    x1 = x1_norm(1,:)';
    y1 = x1_norm(2,:)';
   
    x2 = x2_norm(1,:)';
    y2 = x2_norm(2,:)';

W= [x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(size(x1,1),1)];


% Compute SVD, and f will be the last column of V
[~,~,V]=svd(W,0);

% Compose F 
F=[V(1,end),V(2,end),V(3,end);
   V(4,end),V(5,end),V(6,end);
   V(7,end),V(8,end),V(9,end)];

% Compute again the SVD of F and force F to be of rank 2, so D has a zero in

[U,S,V]=svd(F);

% the his last value

F=U*[S(1,1),0,0;
     0,S(2,2),0;
     0 , 0  , 0]*V';

% Denormalization

F= T2'*F*T1;

end

