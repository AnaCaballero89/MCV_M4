function [Pproj, Xproj] = factorization_method(x1,x2, init) %Pproj, Xproj
    % Outline of the Algoritm
    % 1 Normalize the image coordinates, by applying transformations
    [normX1, T1] = normalise2dpts(x1); % normalises 2D homogeneous points
    [normX2, T2] = normalise2dpts(x2);
    
    lambda = ones (2, size(x1,2));
    
    % Camera1
    % 2. Estimate the fundamental matrices and epipoles
    f1 = fundamental_matrix(x1, x1);
    [~,~, V] = svd(f1);
    e1 = V(:,3)/V(3,3);

    % Camera2
    % 2. Estimate the fundamental matrices and epipoles
    f2 = fundamental_matrix(x2, x1);
    [~,~, V] = svd(f2);
    e2 = V(:,3)/V(3,3);

    for j = 1:size(x1,2)
        % cross = cross product
        %  3. Determine the scale factors
        num = (x1(:,j)'*f1*cross(e1, x1(:,j)));
        denom = (norm(cross(e1,x1(:,j))).^2*lambda(1,j));
        lambda(1,j) = num/denom;
    end
    for j = 1:size(x2,2)
        % cross = cross product
        %  3. Determine the scale factors
        num = x1(:,j)'*f2*cross(e2, x2(:,j));
        denom = norm(cross(e2,x2(:,j))).^2*lambda(1,j);
        lambda(2,j) = num/denom;
    end
   
    
    d= Inf;
    flag = true;
    
    while flag    
        rescale = true;
        i = 0;
        lambdaDiff = Inf;
        %  4. Build the rescaled measurement matrix
        while rescale
            old_lambda_diff = lambdaDiff;
            old_lambda = lambda;
            if mod(i, 2) % normalize rows
                lambda(1,:) = lambda(1,:)./norm(lambda(1,:));
                lambda(2,:) = lambda(2,:)./norm(lambda(2,:));
            else % normalize columns
                for col = 1:size(x1,2)
                    lambda(:,col) = lambda(:,col)/norm(lambda(:,col));
                end
            end           
            % compute euclidan difference with old lambda
            lambdaDiff = (old_lambda - lambda).^2;
            lambdaDiff = sum(lambdaDiff(:));
            lambdaDiff = sqrt(lambdaDiff);
                        
            if ((abs(lambdaDiff - old_lambda_diff)/lambdaDiff) < 0.1)
                rescale = false;
            end
            i = i +1;
        end  % end while rescale
        
        % 5. Balance ,by column-wise and “triplet-of-rows”-wise scalar mutliplication
        M = zeros(3*2, size(x1,2));
        M(1,:) = lambda(1,:).*normX1(1,:);
        M(2,:) = lambda(1,:).*normX1(2,:);
        M(3,:) = lambda(1,:).*normX1(3,:);
        M(4,:) = lambda(2,:).*normX2(1,:);
        M(5,:) = lambda(2,:).*normX2(2,:);
        M(6,:) = lambda(2,:).*normX2(3,:);
        
        % 6. Compute the SVD of the balanced matrix
        [U,D,V] = svd(M);
        
        Pm = U * D(:,1:4);
        Xproj = V(:,1:4)';
        
        dOld= d;
        d = 0;
        
        % 7. From the SVD, recover projective motion and shap
        for i=1:2
            if (i==1)
                 Px = Pm(1:3,:) * Xproj;
                 x = normX1;
            elseif (i==2)
                 Px = Pm(4:6,:) * Xproj;
                 x = normX2;
            end
            for j=1:size(x1,2)
                 d = d + sum((x(:,j) - Px(:,j)).^2);
            end
        end
        
        if ((abs(d - dOld)/d) < 0.1)
            flag = false;
        else % If doesn't converge, update lambdas
            aux = Pm*Xproj;
            lambda(1,:) = aux(3,:);
            lambda(2,:) = aux(6,:);
        end
    end
    
    % 8. Adapt projective motion, to account for the normalization transformations '  of step
    Pproj(1:3,:) = inv(T1)*Pm(1:3,:);
    Pproj(4:6,:) = inv(T2)*Pm(4:6,:);

end