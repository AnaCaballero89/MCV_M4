function [ xnorm, T ] = normalise2( x )
    x(1,1:4)=x(1,1:4)./x(3,1:4);
    x(2,1:4)=x(2,1:4)./x(3,1:4);
    x(3,1:4)=1;
    
    centroid = mean(x(1:2,1:4)')';
    new_x(1,1:4) =  x(1,1:4) - centroid(1);
    new_x(2,1:4) = x(2,1:4) - centroid(2);
    
    dist = sqrt(new_x(1,1:4).^2 + new_x(2,1:4).^2);
    meandist = mean(dist(:));  
    
    s = sqrt(2)/meandist; % scale
    
    T = [s  0   -s*centroid(1)
         0  s   -s*centroid(2)
         0  0   1];
    
     xnorm = T*x;
end