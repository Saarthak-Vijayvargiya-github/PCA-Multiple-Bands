clc;
clear;
close all;

% Specify the bands inside the specBands as 2D Matrices
% NOTE: All the bands should be of same size

specBands = {
    [255, 100; 150, 200];
    [200, 50; 100, 50];
    [50, 150; 200, 255];
    };

% Example for band images
% specBands = {
%     imread("sample.jpg");
%     imread("sample2.jpg");
%     };

bands = cat(3,specBands{1:end});


% Number of bands
dim = size(bands,3);
[M,N] = size(bands(:,:,1));

% Band averages
band_bars = zeros(1,dim);
band_bars(:) = mean(bands(:,:,1:end),[1 2]);

% Band column vectors
band_vecs = zeros(numel(bands(:,:,1)),dim);
band_vecs(:,1:end) = reshape(bands(:,:,1:end),numel(bands(:,:,1)),dim);

% Calculating covariances
covarianceMatrix = zeros(dim,dim);
for i=1:dim
    covarianceMatrix(i,:) = (1/(numel(bands(:,:,i))-1))*sum((band_vecs(:,i)-band_bars(i)).*(band_vecs(:,1:end)-band_bars(1:end)));
end

if(dim<11)
    disp("Covariance Matrix: ");
    disp(covarianceMatrix);
end

% Calculating Eigenvalues and Eigenvectors matrix
eigenVal = eig(covarianceMatrix);
[transFormMat,D] = eig(covarianceMatrix);
if(dim<11)
    disp("Transformation matrix:");
    disp(transFormMat);
end

AllBandMat = band_vecs(:,1:dim) - band_bars(1:dim);
PCA = AllBandMat*transFormMat;

% Getting Principal Components
pCs = [];
PrinComps = cell(2,dim);
for i = 1:dim
    pca_temp = reshape(PCA(:,i),size(bands(:,:,i)));
    pCs = cat(3,pCs,pca_temp);
    PrinComps{1,i} = round(eigenVal(i),3);
    PrinComps{2,i} = pCs(:,:,i);
end

for i = 1:dim
    fprintf("PC%d: eigenVal = %f\n",i,eigenVal(i));
    if(size(bands(:,:,i))<11), disp(pCs(:,:,i)); end
end

% PCA Averages [comes as zero always]
% pca_bars = zeros(1,dim);
% for i = 1:dim
%     pca_bars(i) = mean2(prinComponents(:,:,i));
% end

% Calculating covariances of principal components
Cov_PCA = zeros(dim,dim);
for i=1:dim
    Cov_PCA(i,:) = round((1/(numel(pCs(:,:,i))-1))*sum((PCA(:,i)).*(PCA(:,1:end))),3);
end

% Checking if the principal components are orthogonal (Optional)
if(dim < 11)
    disp("Covariance of PCA:");
    disp(Cov_PCA);
end

% Example 2
% specBands = {
%     [3 3 5 6
%      3 4 4 5
%      4 5 5 6
%      4 5 5 6];
% 
%     [3 2 3 4
%      1 5 3 6
%      4 5 3 6
%      2 4 4 5];
%      
%     [4 2 3 4
%      1 4 2 4
%      4 3 3 5
%      2 3 5 5];
%      };