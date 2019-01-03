function distance = l2_dist(vec) 
% --> needed for randsac_homography 
% --> compute the symmetric geometric error
    distance = sqrt(sum(vec.^2));