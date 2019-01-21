function [disparity] = stereo_computation(left_im, right_im, minDisp, maxDisp, ws, mc)
% input : 
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
% - weight

    disparity = zeros(size(left_im));    
    [h,w] = size(left_im);
    
    pad = floor(ws/2);
    
    left_im = padarray(left_im, [pad pad]);
    right_im = padarray(right_im, [pad pad]);
    
    for i=1+pad:(h + pad)
        for j=1+pad:(w + pad)
            left_patch = double(left_im((i - pad):(i + pad), (j - pad):(j + pad)));
            % size(leftPatch)
            minCol = max((1 + pad), (j - maxDisp));
            maxCol = min((w + pad), (j + maxDisp));
            
            % define window weight
            if strcmp(mc, 'BW')
                
            else    
                weight = zeros(size(left_patch));
                weight(:) = 1/(prod(size(left_patch)));

            if strcmp(mc, 'SSD')
                bestSSD = Inf; 
                for k = minCol:maxCol
                    right_patch = double(right_im(i-pad:i+pad, k-pad:k+pad));
                    ssd = sum(sum( weight.*(left_patch-right_patch).^2 ));
                    if ssd < bestSSD
                        bestSSD = ssd;
                        best = k;
                    end
                end
            elseif strcmp(mc, 'NCC')
                bestNCC = -Inf;                
                for k = minCol:maxCol
                    right_patch = double(right_im(i-pad:i+pad, k-pad:k+pad));
                    sumLeft = sum(left_patch(:).*weight(:));
                    sumRight = sum(right_patch(:).*weight(:));
                    
                    sigmaLeft = sqrt(sum( weight(:).* (left_patch(:) - sumLeft).^2 ));
                    sigmaRight = sqrt(sum( weight(:).* (right_patch(:) - sumRight).^2 ));
                    
                    ncc = sum( weight(:).*(left_patch(:)-sumLeft).*(right_patch(:)-sumRight) )/(sigmaLeft*sigmaRight);
                    if ncc > bestNCC
                        bestNCC = ncc;
                        best = k;
                    end
                end
                
             elseif strcmp(mc, 'BW') 
                 bestBW = Inf; 
                 g_c=14;%gamma c
                 g_p=ws/2;%gamma p,field of view of the human visual system
                 T=20;
                 for k = minCol:maxCol
                    E=0;
                    right_patch = double(right_im(i-pad:i+pad, k-pad:k+pad));
                    for i=1:ws
                        for j=1:ws
                            %Euclidian distance spatial and in the color
                            cp_dist1=abs(left_patch(i,j)-left_patch(pad,pad));
                            cp_dist2=abs(right_patch(i,j)-right_patch(pad,pad));
                            gp_dist=sqrt((i-tab)^2+(j-tab)^2);
                            %Weights 
                            wpq1=exp(-(cp_dist1./g_c+gp_dist./g_p));
                            wpq2=exp(-(cp_dist2./g_c+gp_dist./g_p));
                            %absolute difference AD
                            e=min(abs(left_patch(i,j)-right_patch(i,j)),T); 
                     
                            %dissimilarity BW
                            E=E+(sum(wpq1.*wpq2.*e)/sum(wpq1.*wpq2));
                        end
                    end
                     
                     
                    if E < bestBW
                        bestBW = E;
                        best = E;
                    end
                 
                 end
            
            end 
            end
            
            disparity(i-pad, j-pad) = abs(j-best);
        end        
    end
    a=1;
end
