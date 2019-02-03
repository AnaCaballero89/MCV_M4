function [Pproj, Xproj] = factorization_method(x1,x2,in) 
    % Outline of the Algoritm
    % Normalize the image coordinates, by applying transformations
    % --- cam1 ---
    [q1, T1] = normalise2dpts(x1); 
    % --- cam2 ---
    [q2, T2] = normalise2dpts(x2); 
    
    
    lamda = ones (2, size(x1,2));


    % Estimate the fundamental matrices and epipoles for each camera
    % --- cam1 ---
    % fundamental matrix
    f1 = fundamental_matrix(x1, x1);
    % epipoles
    % https://books.google.es/books?id=cQ9uCQAAQBAJ&pg=PA563&lpg=PA563&dq=svd+epipoles&source=bl&ots=GtE6K9aQxL&sig=ACfU3U0dNHtx1ntjntYQcYpyC9b-tIYGiQ&hl=ca&sa=X&ved=2ahUKEwicwdO67pPgAhWtgM4BHUhdCMsQ6AEwAnoECAQQAQ#v=onepage&q=svd%20epipoles&f=false
    [~,~, v] = svd(f1);
    e1 = v(:,3)/v(3,3);

    % --- cam2 ---
    % fundamental matrix
    f2 = fundamental_matrix(x2, x1);
    % epipoles
    [~,~, v] = svd(f2);
    e2 = v(:,3)/v(3,3);

    % Determine the scale factors using eq 3 (lambda_ip)
    % starting from some arbitrary initial value such as lambda_1p=1

    % a = (e_ij ^ q_ip)*(F_ij q_jp)
    % b = ||e_ij ^q_ip||^2
    % lambda_ip = a/b * lambda_jp
    %% --- cam1 ---
    for p = 1:size(x1,2)   
        q1_ip = x1(:,p);
        q2_ip = x2(:,p);
        e1_ij = e1;
        f1_ij = f1;

        a = (q2_ip'*f1_ij)*cross(e1_ij, q1_ip)  ;
        b = cross(e1_ij,x1(:,p));
        lamda(1,p) = a/norm(b).^2 * lamda(1,p);
    end
    % --- cam2 ---
    for p = 1:size(x2,2)
        q1_ip = x1(:,p);
        q2_ip = x2(:,p);
        e2_ij = e2;
        f2_ij = f2;

        a = (q1_ip'*f2_ij)*cross(e2_ij, q2_ip);
        b = cross(e2_ij,q2_ip);
        lamda(2,p) = a/norm(b).^2 * lamda(1,p);
    end

 
    
    
    % Balance W ,by column-wise and ?triplet-of-rows?-wise scalar mutliplication,
     d= Inf;
     converge = false;
    
     while (converge == false)
 
        %init vars
        repeat = true;
        i = 0;
        lambdaDiff = Inf;
        
        while repeat     
            %auxLamda and auxLamdaDiff are the old values    
            auxLambda = lamda;
            auxLambdaDiff = lambdaDiff;
            
            if mod(i, 2) % Rescale each triplet of rows, pg 5 of the paper, (2)
                lamda(1,:) = lamda(1,:)./norm(lamda(1,:));
                lamda(2,:) = lamda(2,:)./norm(lamda(2,:));
            else %  Rescale each colum, pg 5 of the paper, (1)
                for col = 1:size(x1,2)
                    lamda(:,col) = lamda(:,col)./norm(lamda(:,col));
                end
            end           
            
            % If the entries of W changed significantly, repeat pg 5 of the paper, (3)
            % We consider a significant change from 10%, 
            % 2% is not enought, 5% it's very tight, 15% it's too much
            % Euclidan difference, lambda and old lambda
            lambdaDiff = sqrt(sum((auxLambda - lamda).^2));
            if ((abs(lambdaDiff - auxLambdaDiff)/lambdaDiff) < 0.1)
                repeat = false;
            end
            i = i +1;
        end  % end while rescale

        if(in==1)
          lamda = ones(2,size(x1,2)); 
        end  

        % Build the rescaled measurement matrix W
        % pg 5 of the paper, (3.2)
        W = zeros(3*2, size(x1,2));
        % --- cam1 ---
        W(1,:) = lamda(1,:).*q1(1,:);
        W(2,:) = lamda(1,:).*q1(2,:);
        W(3,:) = lamda(1,:).*q1(3,:);
        % --- cam2 ---
        W(4,:) = lamda(2,:).*q2(1,:);
        W(5,:) = lamda(2,:).*q2(2,:);
        W(6,:) = lamda(2,:).*q2(3,:);
        
        % Check: Notice that with the correct projective depths lambda, the 3mxn measurement matrix W has rank at most rank 4
        size(W);
    
        % Compute the SVD of the balanced matrix W
        [u,s,v] = svd(W);
        
        % From the SVD, recover projective motion and shape
        Pmotion = u*s(:,1:4);
        Xproj = v(:,1:4)';
        
        d_old = d;%d_old is the convergence value of the prev iteration
        d = 0;
        
        % --- cam1 ---
        Pshape = Pmotion(1:3,:) * Xproj;
        x = q2;
        for p=1:size(x1,2)
             d = d + sum((x(:,p) - Pshape(:,p)).^2);
        end
       
        % --- cam2 ---
        Pshape = Pmotion(4:6,:) * Xproj;
        x = q2;
        for p=1:size(x1,2)
             d = d + sum((x(:,p) - Pshape(:,p)).^2);
        end
       
        % we consider a 10% as before for convergence
        if ((abs(d - d_old)/d) < 0.1)
            converge = true;
        else 
            aux = Pmotion*Xproj;
            lamda(1,:) = aux(3,:);
            lamda(2,:) = aux(6,:);
        end
    end
    
    % 8. Adapt projective motion, to account for the norm transf T of step 1
    Pproj(1:3,:) = T1\Pmotion(1:3,:);  % --- cam1 ---
    Pproj(4:6,:) = T2\Pmotion(4:6,:);  % --- cam2 ---

end