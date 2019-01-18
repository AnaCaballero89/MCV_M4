function [X] = triangulate(x1, x2, P1, P2, imsize)
    % Euclidian to Homogeneous coords
    x1=[xe; ones(1, size(x1, 2))];
    x2=[xe; ones(1, size(x2, 2))];
    
    % build matrix H, slide 10 lecture 7
    %  H = [2/nx    0       -1;
    %       0       2/ny    -1;
    %       0       0       1]
      
    H = [2/imsize(1)    0           -1; 
        0               2/imsize(2) -1; 
        0               0           1];
    
    % Transform x, x', P, P' using H  (non-parallel views)
    x1 = euclid(H*x1);  % x
    x2 = euclid(H*x2);  % x'
    P1 = H*P1;          % P
    P2 = H*P2;          % P'
    
    x1p1_3t = x1(1)*P1(3,:);
    p1_1t = P1(1,:);
    
    x2p2_3t = x2(1)*P2(3,:);
    p2_1t = P2(1,:);
    
    y1p1_3t = x1(2)*P1(3,:);
    p1_2t = P1(2,:);
      
    y2p2_3t=x2(2)*P2(3,:);
    p2_2t=P2(2,:);
    
    
    A=[ x1p1_3t - p1_1t;
        y1p1_3t - p1_2t;
        x2p2_3t - p2_1t;
        y2p2_3t - p2_2t ];
    
    
    [~,~, V] = svd(A);
    
    X= V(:,size(V,1));
    X = X ./X(end); % Set 4th coordinate to 1
end