function [p_new,T] = normalise_pts(p)

p(1,:) = p(1,:)./p(3,:);
p(2,:) = p(2,:)./p(3,:);
p(3,:) = 1;

c = mean(p(1:2,:)')';

p_new(1,:)= p(1,:)-c(1);
p_new(2,:)= p(2,:)-c(2);

scale= sqrt(2)/mean(sqrt(p_new(1,:).^2 + p_new(2,:).^2));

T= [scale   0  -scale*c(1)
      0  scale -scale*c(2)
      0     0       1     ];

p_new=T*p;