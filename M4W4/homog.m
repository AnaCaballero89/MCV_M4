function [ xh ] = homog( xe )
    xh = [xe; ones(1, size(xe, 2))];
end