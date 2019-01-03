function  [H] = homography2d(x1, x2)
    [x1_norm, T1] = normalise2(x1);
    [x2_norm, T2] = normalise2(x2);
    x2 = x2_norm(1,:);
    y2 = x2_norm(2,:);
    w2 = x2_norm(3,:);
    
    A = [];
    for i=1:size(x1_norm,2)
        A = [A; 
            zeros(3,1)'          -w2(i)*x1_norm(:,i)'    y2(i)*x1_norm(:,i)'; 
            w2(i)*x1_norm(:,i)'   zeros(3,1)'            -x2(i)*x1_norm(:,i)'];
    end
    
    [~,~,V] = svd(A);
    H = reshape(V(:,9),3,3)';
    H = inv(T2) * H * T1;
end