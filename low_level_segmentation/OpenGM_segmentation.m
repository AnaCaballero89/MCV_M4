clear all;
close all;
clc;

im_name='3_12_s.bmp';

% TODO: Update library path
%Add  library paths

basedir='~/OpenGM_examples/Matlab/';

addpath(basedir);
addpath([basedir,'model/']);
addpath([basedir,'model/functions/']);
addpath([basedir,'../']);

addpath('/usr/local/matlab/mex/');


%Set model parameters
%cluster color
K=10; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[0.0 100]; % Potts Model

%Load images
im = imread(im_name);


%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%
% im = RGB2Lab(im);



%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term
% data_term = -log P( color at pixel 'x' | Cluster color 'c' )  


data_term=[];




%Building 4-grid
%Build OpemGM Model for 4-connected segmentation
disp('create OpenGM model');
% Create OpenGM data
%
gm = CreateGridOpenGMModel(data_term,smooth_term);


if gm ~= -1

    % color clustering
    [~,c] = min(reshape(data_term,[size(im,1)*size(im,2) K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different OpenGM inference algorithmss
    
    inferenceResults = opengm('model', gm, 'a', 'TRBP', 'p', 1, 'v', 'maxIt', 2, 'bound', 0.01,'damping',0.5);
    im_trbp= reshape(mu_color(inferenceResults.states+1,:),size(im));
    
    inferenceResults = opengm('model',gm,'a','ALPHABETASWAP');
    im_ab= reshape(mu_color(inferenceResults.states+1,:),size(im));
    
    % TODO: apply other inference algorithms and compare their performance
    %
    % - LoopyBP
    % - Linear Programing Relaxation
    
    
    figure
    subplot(2,2,1),imshow(Lab2RGB(im));xlabel('Original');
    subplot(2,2,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
    subplot(2,2,3),imshow(Lab2RGB(im_trbp),[]);xlabel('TRBP');
    subplot(2,2,4),imshow(Lab2RGB(im_ab),[]);xlabel('Graph Cut \alpha-\beta swap');
    
else
    
    error('You have to implement the CreateGridOpenGMModel.m function');
    
end