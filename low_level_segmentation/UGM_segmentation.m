clear all;
close all;
clc;

im_name='2_1_s.bmp';

%Set model parameters
%cluster color
K=4; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[0.0 2]; % Potts Model

%Load images
im = imread(im_name);

NumFils = size(im,1);
NumCols = size(im,2);


%Preparing data for GMM fiting
%
% define the unary energy term: data_term
% nodePot = P( color at pixel 'x' | Cluster color 'c' )
im=double(im);
x=reshape(im,[size(im,1)*size(im,2) size(im,3)]);
gmm_color = gmdistribution.fit(x ,K);
mu_color=gmm_color.mu;

data_term = gmm_color.posterior(x);
nodePot = data_term;


%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
[edgePot,edgeStruct] = CreateGridUGMModel(NumFils, NumCols, K ,smooth_term);

if ~isempty(edgePot)

    % color clustering
    [~,c] = max(reshape(data_term,[NumFils*NumCols K]),[],2);
    im_c = reshape(mu_color(c,:),size(im))/255;
    
    % Call different UGM inference algorithms
    display('Loopy Belief Propagation'); tic;
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);toc;
    [~, a] = max(nodeBelLBP,[],2);
    im_lbp = reshape(mu_color(a,:), size(im))/255;
    
    % Max-sum
    display('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    im_bp = reshape(mu_color(decodeLBP,:),size(im))/255;
    toc;
    
    
    % apply other inference algorithms and compare their performance
    % MaxOfMarginals
    disp('Max of Marginals'); tic;
    decodeMaxOfMarg = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_LBP);
    im_MaxOfMarg = reshape(mu_color(decodeMaxOfMarg,:), size(im))/255;
    toc;
    
    % ICM
    disp('ICM'); tic;
    decodeICM = UGM_Decode_ICM(nodePot,edgePot,edgeStruct);
    im_icm = reshape(mu_color(decodeICM,:), size(im))/255;
    toc;
    

    figure
    subplot(2,3,1),imshow(im/255);xlabel('Original');
    subplot(2,3,2),imshow(im_c);xlabel('Clustering without GM');
    subplot(2,3,3),imshow(im_bp);xlabel('Max-Sum');
    subplot(2,3,4),imshow(im_lbp);xlabel('Loopy Belief Propagation');
    
    subplot(2,3,5),imshow(im_MaxOfMarg);xlabel('Max of Marginals');
    subplot(2,3,6),imshow(im_icm);xlabel('ICM');
    
else
   
    error('You have to implement the CreateGridUGMModel.m function');

end