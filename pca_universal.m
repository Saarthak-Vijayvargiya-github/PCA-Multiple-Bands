clc;
clear;
close all;

% Specify the bands inside the specBands as 2D Matrices
% NOTE: All the bands should be of same size
specBands = {
    [255, 100; 150, 200];
    [200, 50; 100, 50];
    [50, 150; 200, 255]
    };

% bands = cat(3,specBands{1:end});

bands = double(imread("sample.jpg"));

% Number of bands
dim = size(bands,3);

% Band averages
band_bars = zeros(1,dim);
for i = 1:dim
    band_bars(i) = mean2(bands(:,:,i));
end

% Band column vectors
band_vecs = zeros(numel(bands(:,:,1)),dim);
for i = 1:dim
    band_vecs(:,i) = reshape(bands(:,:,i),numel(bands(:,:,i)),1);
end

% Calculating covariances
covarianceMatrix = zeros(dim,dim);
for i=1:dim
    for j=1:dim
        covarianceMatrix(i,j) = (1/(numel(bands(:,:,i))-1))*sum((band_vecs(:,i)-band_bars(i)).*(band_vecs(:,j)-band_bars(j)));
    end
end

if(dim<11)
    disp("Covariance Matrix: ");
    disp(covarianceMatrix);
end

eigenVal = eig(covarianceMatrix);
[transFormMat,D] = eig(covarianceMatrix);
if(dim<11)
    disp("Transformation matrix:");
    disp(transFormMat);
end

RGBMat = band_vecs(:,1:dim) - band_bars(1:dim);
PCA = RGBMat*transFormMat;

% Getting Principal Components
prinComp = [];
for i = 1:dim
    pca_temp = reshape(PCA(:,i),size(bands(:,:,i)));
    prinComp = cat(3,prinComp,pca_temp);
end

for i = 1:dim
    fprintf("PC%d: eig = %f\n",i,eigenVal(i));
    if(size(bands(:,:,i))<11), disp(prinComp(:,:,i)); end
end

% PCA Averages come as zero always
% pca_bars = zeros(1,dim);
% for i = 1:dim
%     pca_bars(i) = mean2(prinComponents(:,:,i));
% end

% Calculating covariances of principal components
Cov_PCA = zeros(dim,dim);
for i=1:dim
    for j=1:dim
        Cov_PCA(i,j) = round((1/(numel(prinComp(:,:,i))-1))*sum((PCA(:,i)).*(PCA(:,j))),3);
    end
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