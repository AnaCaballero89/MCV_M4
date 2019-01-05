function [ Y ] = gs_errfunction( P0, Xobs)
    %gs_errfunction 

    %Get H
    H = P0(1:9);
    H = reshape(H,[3 3]);
    
    %Get x of P0
    x = P0(10:end);
    x = reshape(x, [2,size(x,1)/2]);%<----
 
    %homogeneous coordinates
    x = [x; ones(1,size(x,2))];
    xp = H*x;
    
    %Get x1 and x2 observed 
    sizex = size(Xobs,1)/2;
    x1 = Xobs(1:sizex);
    x1 = reshape(x1, [2,size(x1,1)/2]); %<----
    x2 = Xobs(sizex+1:end); 
    x2 = reshape(x2, [2,size(x2,1)/2]); %<----

% homogeneous coordinates
%     xrsize=length(x)/2;
%     xr1=x(1:xrsize);
%     xr2=x(xrsize+1:end);
%     %%xre=[xr1,xr2,ones(xrsize,1)];
%     xrproj=H*xre';
% homogeneous coordinates
%     x=[x1(1:length(x1)/2),x1((length(x1)/2)+1:end),ones(length(x1)/2,1)];
%     xp=[x2(1:length(x2)/2),x2((length(x2)/2)+1:end),ones(length(x2)/2,1)];

    %Error
    vec1 =  x1-euclid(x);
    vec2 = x2-euclid(xp);
    
    Y= (sum(vec1.^2))+(sum(vec2.^2)); %<----
    %%Y=sum((x-xre).^2)+sum((xp-xrproj').^2);

end
