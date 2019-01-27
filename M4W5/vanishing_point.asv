function [v1] = vanishing_point(xo1, xf1, xo2, xf2)
    l1 = cross(xo1,xf1);    % cross product
    l1 = l1/l1(end);        % normalize

    l2 = cross(xo2,xf2);
    l2 = l2/l2(end);

    v1 = cross(l1,l2); 
end