%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
% im1rgb = imread('Data/im1.jpg');  % im2rgb = imread('Data/im2.jpg');
% im3rgb = imread('Data/im3.jpg');  % im4rgb = imread('Data/im4.jpg');

% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.

i1 = imread('Data/im1.jpg');
i2 = imread('Data/im2.jpg');
[matches, p1, p2] = matchOpcional(i1, i2);    % 1: [ f1, f2] ? Match(I1,I2);
[fD1, fS1, fDk, fSk] = Classify_Dyn_Stat(p1, p2);                         % 2: [ fD1 , fS1 , fD2 , fS2 ] = Classify_Dyn_Stat_Ref( f1, f2)
N=4;
for k = N-2                                                                 % 3: for each Ik and k = 3 to N do
    [matches, p1, pk] = matchOpcional(imread('Data/im1.jpg'), imread(strcat('Data/im',int2str(K+2),'.jpg'))); % 4: [ f1, fk ] ? Match(I1,Ik );
    [fDk, fSk] = Classify_Dyn_Stat2(p1, pk);                      % 5: [ fDk , fSk ] = Classify_Dyn_Stat( fS1 , fD1 , fk )
    [F, inliers] = ransac_fundamental_matrix(fDk, fSk, 2.0);                  % 6: Fk =ComputeFundamentalMat( fS1 , fSk ). 
end                                                                         % 7: end for

for df =...;                                                                % 8: for each dynamic feature df ? f1 do 
    l= p1 * p2;                                                             % 9: l = p1 x p2 (l =img line)
    for pk =...;                                                            % 10: for each pk (tk ) ? Si do
        lk= ...*;                                                         % 11: lk=  Fk * pk (l =img line)
        pk = l*lk;                                                          % 12: pk(tk ) = l x lk (pk is the intersection point)
        a = ComputeAlpha(pl p2, pk );                                       % 13: ?k ? ComputeAlpha(pl p2, p(tk ))
    end                                                                     % 14: end for
    o = sort(a);                                                            % 15: ?i ? sort(?k)
end                                                                         % 16: end for


% Thus, the static features are easily
% detected by thresholding the Euclidean distance between
% matched features, and the dynamic features are the remaining ones. 

function [fD1, fSk, fD2, fSk] = Classify_Dyn_Stat(f1, f2)    
end 

function [fSk, fSk] = Classify_Dyn_Stat2(fS1, fD1)
end 