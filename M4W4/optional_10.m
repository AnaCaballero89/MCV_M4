%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  10. New view synthesis

% In this task you are asked to implement part of the new view synthesis method
% presented in the following paper:
% S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.

% You will use a pair of rectified stereo images (no need for prewarping
% and postwarping stages) and their corresponding ground truth disparities
% (folder "new_view").
% Remember to take into account occlusions as explained in the lab session.
% Once done you can apply the code to the another pair of rectified images 
% provided in the material and use the estimated disparities with previous methods.

close all; clear all;

%% Read images
I0 = imread('Data/new_view/im0.png');
I1 = imread('Data/new_view/im1.png');

% http://www.cs.virginia.edu/~cab6fh/pfm.html
disp0 = parsePfm("Data/new_view/disp0.pfm"); 
disp1 = parsePfm("Data/new_view/disp1.pfm");

% Original img
figure; subplot(1,2,1); imshow(I0);
        subplot(1,2,2); imshow(I1);

% Disparities
figure; subplot(1,2,1); imshow(uint8(disp0));
        subplot(1,2,2); imshow(uint8(disp1));        

% Find leftmost and rightmost column positions to generate
% the new empty "morphed" image to put the results in the end
minCol = 0;
maxCol = 0;
for row=1:size(I0,1)
    for col = 1:size(I0,2)
       if (col-disp1(row,col) < minCol) && disp1(row,col) ~= inf
           minCol = col-disp1(row,col);
       elseif (col-disp1(row,col) > maxCol) && disp1(row,col) ~= inf
           maxCol = col-disp1(row,col);
       end
    end 
end

shiftLeft = abs(round(minCol));
shiftRight =  size(I0,2)-round(maxCol);

e = 1;  % Left-right consistency check
s = 0.5;
%for s = 0:0.5:1 % Try different s val
    newView = zeros(size(I0,1), size(I0,2)- shiftLeft -shiftRight, 3); 
    % Obtain the new image
    for row = 1:size(I0,1) 
        for col = shiftLeft:(size(I0,2)-shiftRight)
            % displacement
            if (disp1(row,col) ~= inf) && round(col-disp1(row,col)) >= 1% Avoid idx <= 0
                %if abs(disp1(row, col)-disp0(row, round(col-disp1(row,col)))) < e % Take care of Occlusions 
                    p = (1-s) * [row,col] + s*[row, col-disp1(row,col)];
                    p = round(p);
                    val = (1-s)*I1(row,col,:)+s*I0(row, col-round(disp1(row,col)),:);
                    newView(p(1),p(2),:) = val;
                %end
            end
        end 
    end
    figure; imshow(uint8(newView));
%end