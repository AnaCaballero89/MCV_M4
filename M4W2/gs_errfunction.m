function [ Y ] = gs_errfunction( P0, Xobs)
%gs_errfunction 

%Get H
H = P0(1:9);
H=reshape(H,[3 3]);
%Get x of P0
xr=P0(10:end);
%homogeneous coordinates

xrsize=length(xr)/2;
xr1=xr(1:xrsize);
xr2=xr(xrsize+1:end);
xre=[xr1,xr2,ones(xrsize,1)];
xrproj=H*xre';

%Get x and xp observed 
sizex=size(Xobs,1)/2;
x1=Xobs(1:sizex);
x2=Xobs(sizex+1:end);
%homogeneous coordinates
x=[x1(1:length(x1)/2),x1((length(x1)/2)+1:end),ones(length(x1)/2,1)];
xp=[x2(1:length(x2)/2),x2((length(x2)/2)+1:end),ones(length(x2)/2,1)];

%Error
Y=sum((x-xre).^2)+sum((xp-xrproj').^2);

end

