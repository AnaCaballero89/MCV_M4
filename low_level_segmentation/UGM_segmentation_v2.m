clear all;
close all;
clc;

im_name='3_12_s.bmp';

% TODO: Update library path = done
% Add  library paths
basedir='~/Desenvolupament/UGM/';
addpath(basedir);



%Set model parameters
%cluster color
K=4; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[0.0 2]; % Potts Model

%Load images
im = imread(im_name);


NumFils = size(im,1);
NumCols = size(im,2);

%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%
im = RGB2Lab(im);



%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term = nodePot
im=double(im);
x=reshape(im,[size(im,1)*size(im,2),size(im,3)]);
gmm_color = gmdistribution.fit(x,K);
mu_color=gmm_color.mu;
nodePot=gmm_color.posterior(x);

%nodePot = P( color at pixel 'x' | Cluster color 'c' )  


%nodePot=[];



%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
[edgePot,edgeStruct] = CreateGridUGMModel(NumFils, NumCols, K ,smooth_term);


if ~isempty(edgePot)

    % color clustering
    [~,c] = min(reshape(data_term,[NumFils*NumCols K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different UGM inference algorithms
    display('Loopy Belief Propagation'); tic;
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);toc;
    im_lbp = max(nodeBelLBP,[],2);
    
    % Max-sum
    display('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    im_bp= reshape(mu_color(decodeLBP,:),size(im));
    toc;
    
    
    % TODO: apply other inference algorithms and compare their performance
    %
    % - Graph Cut
    display('Graph Cut'); tic;
    decodeGP = UGM_Decode_GraphCut(nodePot, edgePot, edgeStruct);
    im_gp= reshape(mu_color(decodeGP,:),size(im));
    toc;
    % - Linear Programing Relaxation
    display('Linear Programing Relaxation'); tic;
    decodeLPR = UGM_Decode_LinProg(nodePot, edgePot, edgeStruct);
    im_lpr= reshape(mu_color(decodeLPR,:),size(im));
    toc;
    
    figure
    subplot(2,3,1),imshow(Lab2RGB(im));xlabel('Original');
    subplot(2,3,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
    subplot(2,3,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
    subplot(2,3,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');

    subplot(2,3,5),imshow(Lab2RGB(im_gp),[]);xlabel('Graph Cut');
    subplot(2,3,5),imshow(Lab2RGB(im_lpr),[]);xlabel('Linear Programing Relaxation');

else
   
    error('You have to implement the CreateGridUGMModel.m function');

end