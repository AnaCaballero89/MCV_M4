function [ F ] = fundamental_matrix(x1_test, x2_test)
% Computation of the Fundamental Matrix with 8 points or more.
%   Detailed explanation goes here
numb=size(x1_test,2); 
%normalization of the points by the function provided 

[x1_norm, ~] = normalise2dpts(x1_test);
[x2_norm, ~] = normalise2dpts(x2_test);


 bu=x1_norm(1,:);
%Compute linear combination of the matches such as W*f=0
W=[x1_norm(1,:)'.* x2_norm(1,:)', x1_norm(2,:)'.* x2_norm(1,:)',...
   x2_norm(1,:)', x1_norm(1,:)'.* x2_norm(2,:)',x2_norm(2,:)',... 
   x1_norm(1,:)', x1_norm(2,:)',ones(1,numb)'];


%Compute SVD, and f will be the last column of V
[U,S,V]=svd(W);
%Compose F 
F=[V(1,end),V(2,end),V(3,end);
   V(2,end),V(4,end),V(5,end);
   V(3,end),V(5,end),V(6,end)];
%Compute again the SVD of F and force F to be of rank 2, so D has a zero in

[U,S,V]=svd(F);

%the his last value

F=U*[S(1,1),0,0;
     0,S(2,2),0;
     0 , 0  , 0]*V';
end

