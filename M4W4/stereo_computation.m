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
                 bestBW =Inf; 
                 g_c=5;%gamma c, given by the paper
                 g_p=ws/2;%gamma p, radius of the window (field of view of the human visual system)
                 T=40;% Truncation values,constant given by the paper
                 center=ceil(ws/2);
                 for k = minCol:maxCol
                    numerador=0;
                    denonminador=0;
                    
                    right_patch = double(right_im(i-pad:i+pad, k-pad:k+pad));
                    for p=1:ws      
                        for q=1:ws
                            %Euclidian distance spatial and in color
                            cp_dist1=(1/3)*abs(double(left_patch(p,q)-left_patch(center,center)));
                            cp_dist2=(1/3)*abs(double(right_patch(p,q)-right_patch(center,center)));
                            gp_dist=sqrt((p-center)^2+(q-center)^2);
                            %Weights 
                            wpq1=exp(-(cp_dist1/g_c+gp_dist./g_p));
                            wpq2=exp(-(cp_dist2/g_c+gp_dist./g_p));
                            %absolute difference AD
                            e=min(abs(double(left_patch(p,q)-right_patch(p,q))),T); 
                     
                            %dissimilarity BW 
                            numerador=numerador+(wpq1*wpq2*e);
                            denonminador=denonminador+(wpq1*wpq2);
                        end
                    end
                    E=numerador/denonminador; 
                    if E < bestBW
                        bestBW = E;
                        best = k;
                    end
                 
                 end
            
              end 
            
            disparity(i-pad, j-pad) = abs(j-best);
        end        
    end

