function [ phi ] = G5_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni )
%Implementation of the Chan-Vese segmentation following the explicit
%gradient descent in the paper of Pascal Getreur "Chan-Vese Segmentation".
%It is the equation 19 from that paper

%I     : Gray color image to segment
%phi_0 : Initial phi
%mu    : mu lenght parameter (regularizer term)
%nu    : nu area parameter (regularizer term)
%eta   : epsilon for the total variation regularization
%lambda1, lambda2: data fidelity parameters
%tol   : tolerance for the sopping criterium
% epHeaviside: epsilon for the regularized heaviside. 
% dt     : time step
%iterMax : MAximum number of iterations
%reIni   : Iterations for reinitialization. 0 means no reinitializacion

[ni,nj]=size(I);
hi=1;
hj=1;


phi=phi_0;
dif=inf;
nIter=0;
while dif>tol && nIter<iterMax
    
    phi_old=phi;
    nIter=nIter+1;        
    
    
    %Fixed phi, Minimization w.r.t c1 and c2 (constant estimation)
    % === eq 10 & 11 IPOL_Getreuer_2012.pdf
    % === the optimal values of c1 and c2 are the region averages
    c1 = sum((phi(:) >= 0) .* I(:))/sum(phi(:) >= 0); %TODO 1: Line to complete
    c2 = sum((phi(:) <  0) .* I(:))/sum(phi(:) <  0); %TODO 2: Line to complete
    
    %Boundary conditions
    % === eq 25
    % === boundary condition is enforced by duplicating pixels near the borders
    phi(1,:)   =  phi_old(2,:); %TODO 3: Line to complete
    phi(end,:) =  phi_old(end-1,:); %TODO 4: Line to complete
    
    phi(:,1)   = phi_old(:,2); %TODO 5: Line to complete
    phi(:,end) = phi_old(:,end-1); %TODO 6: Line to complete

    
    %Regularized Dirac's Delta computation
    delta_phi = G5_diracReg(phi, epHeaviside);   %notice delta_phi=H'(phi)	
    
    %derivatives estimation
    %i direction, forward finite differences
    phi_iFwd  = DiFwd(phi); %TODO 7: Line to complete
    phi_iBwd  = DiBwd(phi); %TODO 8: Line to complete
    
    %j direction, forward finitie differences
    phi_jFwd  = DjFwd(phi); %TODO 9: Line to complete
    phi_jBwd  = DjBwd(phi); %TODO 10: Line to complete
    
    %centered finite diferences
    phi_icent   = (phi_iFwd + phi_iBwd)/2; %TODO 11: Line to complete
    phi_jcent   = (phi_jFwd + phi_jBwd)/2; %TODO 12: Line to complete
    
    %A and B estimation (A y B from the Pascal Getreuer's IPOL paper "Chan
    %Vese segmentation
    % === eq 18
    % ===  curvature is regularized by including a small positive parameter (mu) in the denominators that prevents division by zero
    A = mu./sqrt(eta^2 + phi_jFwd.^2 + phi_icent.^2); %TODO 13: Line to complete
    B = mu./sqrt(eta^2 + phi_jcent.^2 + phi_iFwd.^2); %TODO 14: Line to complete
    
    
    %%Equation 22, for inner points
    %TODO 15: Line to complete
    % === eq 22
    % === phi to be solved by one Gauss–Seidel swee
    numerator =  (phi(2:end-1,2:end-1) + dt*delta_phi(2:end-1,2:end-1) .* ...
        (A(2:end-1,2:end-1) .* phi(3:end,2:end-1) + A(1:end-2,2:end-1) .* phi(1:end-2,2:end-1) + ...
        B(2:end-1,2:end-1) .* phi(2:end-1,3:end) + B(2:end-1,2:end-1) .* phi(2:end-1,1:end-2) + ...
        -nu - lambda1 .* (I(2:end-1,2:end-1) - c1).^2 + lambda2 .* (I(2:end-1,2:end-1) - c2).^2));   
    denumerator1= 1 + dt .* delta_phi(2:end-1,2:end-1);
    denumerator2 = A(2:end-1,2:end-1) + A(1:end-2,2:end-1) + ...
        B(2:end-1,2:end-1) + B(2:end-1,1:end-2);
    phi(2:end-1,2:end-1) = numerator./ (denumerator1 .* denumerator2);
    
    
    %Reinitialization of phi
    if reIni>0 && mod(nIter, reIni)==0
        indGT = phi >= 0;
        indLT = phi < 0;
        
        phi=double(bwdist(indLT) - bwdist(indGT));
        
        %Normalization [-1 1]
        nor = min(abs(min(phi(:))), max(phi(:)));
        phi=phi/nor;
    end
  
    %Diference. This stopping criterium has the problem that phi can
    %change, but not the zero level set, that it really is what we are
    %looking for.
    dif = mean(sum( (phi(:) - phi_old(:)).^2 ));
          
    %Plot the level sets surface
    subplot(1,2,1) 
        %The level set function
        surfc(phi(2:end-1,2:end-1)) 
        %TODO 16: Line to complete 
        hold on
        %The zero level set over the surface
        contour(phi(2:end-1,2:end-1)); %TODO 17: Line to complete
        hold off
        title('Phi Function');
    
    %Plot the curve evolution over the image
    subplot(1,2,2)
        imagesc(I);        
        colormap gray;
        hold on;
        contour(phi(2:end-1,2:end-1)) %TODO 18: Line to complete
        title('Image and zero level set of Phi')

        axis off;
        hold off
    drawnow;
    pause(.0001); 
end