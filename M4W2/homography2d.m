function H = homography2d(x1, x2)

    [x1, T1] = normalise_pts(x1);
    [x2, T2] = normalise_pts(x2);
    
    A = [];
    
    for n = 1:length(x1)
        X = x1(:,n)' ;
        x = x2(1,n);
        y = x2(2,n);
        w = x2(3,n);
        A(3*n-2,:) = [ zeros(1,3)  -w*X   y*X];
        A(3*n-1,:) = [ w*X   zeros(1,3)  -x*X];
        A(3*n  ,:) = [-y*X   x*X   zeros(1,3)];
    end
    
    [~,~,V] = svd(A);

    H = reshape(V(:,9),3,3)';
    
    H = T2\H*T1;
    