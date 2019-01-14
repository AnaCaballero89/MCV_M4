function [matches, p1, p2] = matchOpcional(img1, img2)
    
    im1 = sum(double(img1), 3) / 3 / 255;
    im2 = sum(double(img2), 3) / 3 / 255;
    [points_1, desc_1] = sift(im1, 'Threshold', 0.01);
    [points_2, desc_2] = sift(im2, 'Threshold', 0.01);
    matches = siftmatch(desc_1, desc_2);
    
    % p1 and p2 contain the homogeneous coordinates of the matches
    p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
    p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];
    
    figure;
    plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');
end