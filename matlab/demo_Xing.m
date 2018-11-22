% Demo of using corpca.m for CORPCA
% This function has written based on Programs from Matlab 
%     Copyright (c) 2017, Huynh Van Luong, version 01, Jan. 24, 2017
%     Multimedia Communications and Signal Processing, University of Erlangen-Nuremberg.
%     All rights reserved.
%
%     PUBLICATION: Huynh Van Luong, N. Deligiannis, J. Seiler, S. Forchhammer, and A. Kaup, 
%             "Incorporating Prior Information in Compressive Online Robust Principal Component Analysis," 
%              in e-print, arXiv, Jan. 2017.
%
% 
%% Initialization 
%addpath('./videos');
%addpath('./videos');
addpath(genpath('./'));
fpath = './videos/bootstrap/';
imagefiles = dir('./videos/bootstrap/*.bmp');      
nfiles = length(imagefiles);    % Number of files found

[h0, w0, c] = size(imread([fpath,imagefiles(1).name]));
resizeScale = 0.75;
h = floor(h0*resizeScale);
w = floor(w0*resizeScale);

% s0(1) = 40; % sparse degree
% s0(2) = s0(1) + 15; % sparse constraint of sparse components
nSI = 3; % Size of number of foreground prior information 
m = h*w; % a number of measurements of reduced data  
% sj = max(round(mean(s0)/2),1); % s_j = || x_t - x_t-1||_0
n = h*w; % Dimension of data vectors
% r = 5; % rank of low-rank components
q = 10; % the number of testing vectors
% d = 100; % the number of training vectors



images = zeros([n,nfiles]);    % image size (120,160,3)

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread([fpath,currentfilename]);
   currentimage = imresize(rgb2gray(currentimage), resizeScale);
   images(:,ii) = currentimage(:);
end
trainData = images(:,1:q);
disp('Image Size: ')
disp(size(currentimage));

%% Generating numerical example data for CORPCA
%     trainData: Traning data
%     M= L + S: Testing data
%     L: Testing low-rank components
%     S: Testing sparse components
% [M, trainData, L, S] = dataGeneration (n, q, d, r, s0, sj);
M = images(:,101:100+q);

%% ***********************************************************************************
% Foreground-background separation by CORPCA
% ************************************************************************************
% Initializing background and foreground prior information via an offline rpca
% addpath(genpath('inexact_alm_rpca'));
RPCA_lambda     = 1/sqrt(size(trainData,1));
[B0, Z0, ~] = inexact_alm_rpca(trainData, RPCA_lambda, -1, 20); 

%% Running one CORPCA demo
%--------------------------------------------------------------------------------------------------
%--- Input: observation yt with m measurements, projection matrix Phi, prior information Btm1, Ztm1
%--- Output: recovered foreground xt, background vt, the prior information updates, Bt, Zt
%--------------------------------------------------------------------------------------------------
%{  
rpMat = randn(m, n);
Phi = rpMat; % Input the measurement matrix     
Btm1 = B0; % Input background prior
Ztm1 = Z0(:, end - nSI + 1 : end); % Input foreground prior
fgs = zeros(h,w,q);
bgs = zeros(h,w,q);
tic;
for t = 1 : q
%t = 1; % Testing fame i = 1  
    fprintf('Testing fame %d at a measurement rate %2.2f \n', t, m/n);

    yt = Phi*M (:,t); % Input observation    
    [xt, vt, Zt, Bt] = corpca(yt, Phi, Ztm1, Btm1); % Performing CORPCA for separation
    % update prior information
    Ztm1 = Zt;
    Btm1 = Bt;
    %fg(:,t) = xt;
    %bg(:,t) = vt;
    fgs(:,:,t) = reshape(xt,[h,w]);
    bgs(:,:,t) = reshape(vt,[h,w]);
    
%     fprintf('Recovered foreground error: %4.8f \n', norm(xt - S(:,t),2)/(norm(S(:,t),2)));
%     fprintf('Recovered background error: %4.8f \n\n', norm(vt - L(:,t),2)/(norm(L(:,t),2)));
end
toc;
fgs = reshape(fgs, [h,w,1,q]);
figure;montage(fgs)
%img = reshape(fg,[120,160]);
%}
%% Running one CORPCA-OF demo
%--------------------------------------------------------------------------------------------------
%--- Input: observation yt with m measurements, projection matrix Phi, prior information Btm1, Ztm1
%--- Output: recovered foreground xt, background vt, the prior information updates, Bt, Zt
%--------------------------------------------------------------------------------------------------
  
%rpMatOF = randn(m, n);
rpMatOF = eye(m, n);
PhiOF = rpMatOF; % Input the measurement matrix     
Btm1OF = B0; % Input background prior
Ztm1OF = Z0(:, end - nSI + 1 : end); % Input foreground prior
fgsOF = zeros(h0,w0,q);
bgsOF = zeros(h0,w0,q);

for t = 1 : q
%t = 1; % Testing fame i = 1  
    fprintf('Testing fame %d at a measurement rate %2.2f \n', t, m/n);
    fprintf('Updating proir by optical flow...\n');
    tic;
    % Compute optical flow
    of12 = computeOF(reshape(Ztm1OF(:,end-2+1),[h,w]), reshape(Ztm1OF(:,end-1+1),[h,w]));
    of13 = computeOF(reshape(Ztm1OF(:,end-3+1),[h,w]), reshape(Ztm1OF(:,end-1+1),[h,w]));
    
    % Motion Compensation
    Ztm1OF(:,end-2+1) = linearOFCompensate(reshape(Ztm1OF(:,end-2+1),[h,w]), of12, 1,  true);
    Ztm1OF(:,end-3+1) = linearOFCompensate(reshape(Ztm1OF(:,end-3+1),[h,w]), of13, 0.5,true);

    fprintf('Optimizing...\n');
    yt = PhiOF*M (:,t); % Input observation    
    [xt, vt, Zt, Bt, beta, Wk] = corpca(yt, PhiOF, Ztm1OF, Btm1OF); % Performing CORPCA for separation
    % update prior information
    Ztm1OF = Zt;
    Btm1OF = Bt;
    %fg(:,t) = xt;
    %bg(:,t) = vt;
    fgsOF(:,:,t) = imresize(reshape(xt,[h,w]), 1/resizeScale);
    bgsOF(:,:,t) = imresize(reshape(vt,[h,w]), 1/resizeScale);
    
    %disp('beta:')
    %disp(beta)
    
%     fprintf('Recovered foreground error: %4.8f \n', norm(xt - S(:,t),2)/(norm(S(:,t),2)));
%     fprintf('Recovered background error: %4.8f \n\n', norm(vt - L(:,t),2)/(norm(L(:,t),2)));
    toc;
end
%disp('Beta:')
%disp(beta)


fgsOF = reshape(fgsOF, [h0,w0,1,q]);
figure;montage(uint8(fgsOF))

%% 
v = reshape(M, [h, w, 1, q]);
%figure;montage(uint8(v));
%figure;montage(uint8(fgsOF.*v));

%% Utility functions
function of = computeOF(img1, img2, method)
    if nargin == 2
        method = 'LK';
    end
    if strcmp(method, 'LK')
        ofestmator = opticalFlowLK('NoiseThreshold',0.009);
        estimateFlow(ofestmator, img1);
        of = estimateFlow(ofestmator, img2);
    end
end
