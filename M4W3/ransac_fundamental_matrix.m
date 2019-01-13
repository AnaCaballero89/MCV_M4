function [F, idx_inliers] = ransac_fundamental_matrix (x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);
max_it=1000;

% ransac
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it
    
    points = randomsample(Npoints, 8); % we use the minimum number of points required
    F = fundamental_matrix(x1(:,points), x2(:,points));
    inliers = compute_inliers(F, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute F from all the inliers
F = fundamental_matrix(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(F, x1, x2, th)
    
    % To compute the inliers we use the assumption that two correspondences 
    % x and x' in two images, I and I', must comply that x'* F * x = 0.
    x1 = x1 ./ repmat(x1(end,:), size(x1,1), 1); % normalise x1;
    x2 = x2 ./ repmat(x2(end,:), size(x2,1), 1);% normalise x2;
    x2Fx1 = zeros(1,length(x1));
    
    for k = length(x1)
        x2Fx1(k)= x2(:,k)'*F*x1(:,k);
    end
    
    % Then we compute the epipolar lines for each correspondence where they
    % should be in the same line or should be very close to be inliers
    % depending on the distance or error (sampson distance).
    
    Fx1= F*x1;
    Fx2= F'*x2;
 
    % Instead of compute the symmetric geometric error we use first-order 
    % geometric error (Sampson distance):
    
    num = x2Fx1.^2;
    denom= Fx1(1,:).^2 + Fx1(2,:).^2 + Fx2(1,:).^2 + Fx2(2,:).^2;
    d2 =  num./denom;
    
    idx_inliers = find(abs(d2) < th);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat