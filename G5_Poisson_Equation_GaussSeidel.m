function [u] = G5_Poisson_Equation_GaussSeidel(f, dom2Inp, param)
%this code is not intended to be efficient. 

% ni = h, num rows
% nj = w, num cols

[ni, nj]=size(f);

%We add the ghost boundaries (for the boundary conditions)
% for the img
f_ext = zeros(ni+2, nj+2);
f_ext(2:end-1, 2:end-1) = f;
% for the mask
dom2Inp_ext =zeros(ni+2, nj+2);
dom2Inp_ext(2:end-1, 2:end-1) = dom2Inp;

% Add padding to the driving parameter
b = zeros(ni+2, nj+2);
b(2:end-1, 2:end-1) = param.driving;

% https://www.youtube.com/watch?v=sO4TWTgO4TI
% We can increase the number of iterations as we consider. As bigger it
% gets closer to the perfect result
for numIter= 1:10
    for i=1:ni
        for j=1:nj
            if dom2Inp_ext(i, j) == 1
                upPix = f_ext(i-1, j);
                rightPix = f_ext(i, j+1);
                downPix = f_ext(i+1, j);
                leftPix = f_ext(i, j-1);
                f_ext(i, j) =  ((upPix + rightPix + downPix + leftPix)-b(i, j)) / 4;                
            end
        end
    end
    f_ext(i-100, j-100)
end
    %Eliminate the ghost boundaries
    u = f_ext(2:end-1, 2:end-1);