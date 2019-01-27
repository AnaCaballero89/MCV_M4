function [Pproj, Xproj] = factorization_method(x1,x2, init) %Pproj, Xproj


    [norm_x1, T1] = normalise2dpts(x1); % normalises 2D homogeneous points
    [norm_x2, T2] = normalise2dpts(x2);
    
    lambda = ones (2, size(x1,2));
    
    if isequal(init, 'sturm')
        % Camera1
        f1 = fundamental_matrix(x1, x1);
        [~,~, V] = svd(f1);
        e1 = V(:,3) / V(3,3);
            
        % Camera2
        f2 = fundamental_matrix(x2, x1);
        [~,~, V] = svd(f2);
        e2 = V(:,3) / V(3,3);
            
        for j = 1:size(x1,2)
            % cross = cross product
            num = (x1(:,j)'*f1*cross(e1, x1(:,j)));
            denom = (norm(cross(e1,x1(:,j))).^2*lambda(1,j));
            lambda(1,j) = num/denom;
        end
        for j = 1:size(x2,2)
            num = x1(:,j)'*f2*cross(e2, x2(:,j));
            denom = norm(cross(e2,x2(:,j))).^2*lambda(1,j);
            lambda(2,j) = num/denom;
        end
    end
    
    d= Inf;
    flag = true;
    
    while flag    
        rescale = true;
        counter = 0;
        lambda_diff = Inf;
        
        while rescale
            old_lambda_diff = lambda_diff;
            old_lambda = lambda;
            if mod(counter, 2)
                % normalize rows
                lambda(1,:) = lambda(1,:)./norm(lambda(1,:));
                lambda(2,:) = lambda(2,:)./norm(lambda(2,:));
            else
                % normalize columns
                for col = 1:size(x1,2)
                    lambda(:,col) = lambda(:,col)/norm(lambda(:,col));
                end
            end           
            % compute euclidan difference with old lambda
            lambda_diff = (old_lambda - lambda).^2;
            lambda_diff = sum(lambda_diff(:));
            lambda_diff = sqrt(lambda_diff);
                        
            if ((abs(lambda_diff - old_lambda_diff)/lambda_diff) < 0.1)
                rescale = false;
            end
            counter = counter +1;
        end  % end while rescale
        
        M = zeros(3*2, size(x1,2));
        M(1,:) = lambda(1,:) .* norm_x1(1,:);
        M(2,:) = lambda(1,:) .* norm_x1(2,:);
        M(3,:) = lambda(1,:) .* norm_x1(3,:);
        M(4,:) = lambda(2,:) .* norm_x2(1,:);
        M(5,:) = lambda(2,:) .* norm_x2(2,:);
        M(6,:) = lambda(2,:) .* norm_x2(3,:);
        
        [U,D,V] = svd(M);
        
        Pm = U * D(:,1:4);
        Xproj = V(:,1:4)';
        
        d_old= d;
        d = 0;
       
        for i=1:2
            if i==1
                 Px = Pm(1:3,:) * Xproj;
                 x = norm_x1;
            elseif i==2
                 Px = Pm(4:6,:) * Xproj;
                 x = norm_x2;
            end
            for j=1:size(x1,2)
                 d = d + sum((x(:,j) - Px(:,j)).^2);
            end
        end
        
        if ((abs(d - d_old)/d) < 0.1)
            flag = false;
        else
            % If it has not converged update lambdas
            temp = Pm*Xproj;
            lambda(1,:) = temp(3,:);
            lambda(2,:) = temp(6,:);
        end
    end
    
    Pproj(1:3,:) = inv(T1) * Pm(1:3,:);
    Pproj(4:6,:) = inv(T2) * Pm(4:6,:);

end