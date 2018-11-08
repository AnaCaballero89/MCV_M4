function gm=CreateGridOpenGMModel(im,lambda)
%
%
% im: imput "image" with unary potentials in the 3rd dimension
% lambda: smoothing factor



% parameter
numVariablesN = size(im,1); % number of variables of first dimension
numVariablesM = size(im,2); % number of variables of second dimension
numLabels = size(im,3);     % number of labels for each variable



tic

numVariables = numVariablesN * numVariablesM;

% create model
gm = openGMModel;

% add variables
gm.addVariables(repmat(numLabels, 1, numVariables));

% add unary factor to each variable
DataTerm=reshape(im,[numVariables  numLabels]);

% Each variable is conditioned by the "observed" variable
for n = 1: numVariables
    unaryFunction = openGMExplicitFunction(numLabels, DataTerm(n,:));
    gm.addFunction(unaryFunction);
    gm.addFactors(unaryFunction, n-1);
end


% Pair-wise feature functions
% binary function
pottsFunction = openGMPottsFunction([numLabels, numLabels], lambda(1), lambda(2));
% add functions
gm.addFunction(pottsFunction);



% add binary factors to create grid structure
% horizontal factors
variablesH = 0 : (numVariables - 1);
variablesH(numVariablesN : numVariablesN : numVariables) = [];
variablesH = cat(1, variablesH, variablesH + 1);

% vertical factors
variablesV = 0 : (numVariables - (numVariablesN + 1));
variablesV = cat(1, variablesV, variablesV + numVariablesN);

% concatenate horizontal and vertical factors
variables = cat(2, variablesH, variablesV);

% add factors
gm.addFactors(pottsFunction, variables);
toc;