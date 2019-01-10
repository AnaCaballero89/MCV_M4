function [ F ] = fundamental_matrix(x1_test, x2_test)
% Computatation of the Fundamental Matrix with 8 points or more.
%   Detailed explanation goes here

%Compute linear combination of the matches such as W*f=0
W=[x1_test(1,:)'.* x2_test(1,:)', x1_test(2,:)'.* x2_test(1,:)',...
   x2_test(1,:)', x1_test(1,:)'.* x2_test(2,:)',x2_test(2,:)' ... 
   x1_test(1,:)', x1_test(2,:)',ones(1,size(x1_test(2,:)))];


%Compute SVD, and f will be the last column of V
[U,S,V]=svd(W);
%Compose F 
F=[V(1,end),V(2,end),V(3,end);
   V(2,end),V(4,end),V(5,end);
   V(3,end),V(5,end),V(6,end)];
%Compute again the SVD of F and force F to be of rank 2, so D has a zero in

[U,S,V]=svd(F);

%the his last value

F=U*[S(1,1),S(1,2),S(1,3);
    S(2,1),S(2,2),S(2,3);
    S(3,1),S(3,2),0]*V';
end

