% demo_xing
%% Read and reshape image
% https://www.mathworks.com/help/matlab/ref/videoreader.read.html
% video = read(v)
% video = read(v,index)
% video = read(___,'native')
vRaw = read(VideoReader('./videos/bootstrap.avi')); % 60          80           3        3055
[h, w, c, nFrame] = size(vRaw);
nClip = 500;
vRawClipGray = squeeze(vRaw(:,:,1,1:nClip));
%vRawClipGray = zeros(h,w,nClip);
%for i=1:nClip
%    vRawClipGray(:,:,i) = rgb2gray(vRawClip(:,:,:,i));
%end
v = double(reshape(vRawClipGray,[h*w, nClip]));
%v1 = reshape(v(:,10),[h,w]);
%imshow(v1)
% vRaw = rgb2gray(read('./videos/bootstrap.avi'));
% v = reshape(v,[])

%% CORPCA
% Initializing background and foreground prior information via an offline rpca
qX = 10; % the number of testing vectors
dX = 10; % the number of training vectors
MX = v(:,end-qX+1:end);
trainDataX = v(:,1:dX);
% trainDataX = MX;

addpath(genpath('inexact_alm_rpca'));
RPCA_lambdaX = 1/sqrt(size(trainDataX,1));
[B0X, Z0X, ~] = inexact_alm_rpca(trainDataX, RPCA_lambdaX, -1, 20); 

s0X(1) = 10; % sparse degree
s0X(2) = s0X(1) + 15; % sparse constraint of sparse components
nSI = 3; % Size of number of foreground prior information 
m = h*w; % a number of measurements of reduced data  
sj = max(round(mean(s0X)/2),1); % s_j = || x_t - x_t-1||_0
nX = h*w; % Dimension of data vectors
r = 10; % rank of low-rank components



rpMatX = randn(m, nX);
PhiX = rpMatX; % Input the measurement matrix     
Btm1X = B0X; % Input background prior, "tm1" for "t minus one"?
Ztm1X = Z0X(:, end - nSI + 1 : end); % Input foreground prior
tic

% init optical flow estimator
fgs = zeros(h,w,qX);
bgs = zeros(h,w,qX);
for t = 1 : qX
%t = 1; % Testing fame i = 1      
    fprintf('Testing fame %d at a measurement rate %2.2f \n', t, m/nX);       
    % Performing CORPCA for separation
    ytX = PhiX*MX (:,t); % Input observation
    [xtX, vtX, ZtX, BtX] = corpca(ytX, PhiX, Ztm1X, Btm1X); % Performing CORPCA for separation
    % update prior information
    Ztm1X = ZtX;
    Btm1X = BtX;
    
    fgs(:,:,t) = reshape(xtX,[h,w]);
    bgs(:,:,t) = reshape(vtX,[h,w]);

    %fprintf('Recovered foreground error: %4.8f \n', norm(xtX - S(:,t),2)/(norm(S(:,t),2)));
    %fprintf('Recovered background error: %4.8f \n\n', norm(vtX - L(:,t),2)/(norm(L(:,t),2)));
    %reset(of12estmator);estimateFlow(of12estmator, reshape(Ztm1X(:,end-2+1),[h,w]));
    %reset(of13estmator);estimateFlow(of13estmator, reshape(Ztm1X(:,end-3+1),[h,w]));
end
toc
figure;montage(fgs)

%% CORPCA - OF
% Initializing background and foreground prior information via an offline rpca
qX = 10; % the number of testing vectors
dX = 10; % the number of training vectors
MX = v(:,end-qX+1:end);
% trainDataX = v(:,1:dX);
trainDataX = MX;

addpath(genpath('inexact_alm_rpca'));
RPCA_lambdaX = 1/sqrt(size(trainDataX,1));
[B0X, Z0X, ~] = inexact_alm_rpca(trainDataX, RPCA_lambdaX, -1, 20); 

s0X(1) = 10; % sparse degree
s0X(2) = s0X(1) + 15; % sparse constraint of sparse components
nSI = 3; % Size of number of foreground prior information 
m = h*w; % a number of measurements of reduced data  
sj = max(round(mean(s0X)/2),1); % s_j = || x_t - x_t-1||_0
nX = h*w; % Dimension of data vectors
r = 10; % rank of low-rank components



rpMatX = randn(m, nX);
PhiX = rpMatX; % Input the measurement matrix     
Btm1X = B0X; % Input background prior, "tm1" for "t minus one"?
Ztm1X = Z0X(:, end - nSI + 1 : end); % Input foreground prior
tic

% init optical flow estimator
%of12estmator = opticalFlowLK('NoiseThreshold',0.009);estimateFlow(of12estmator, reshape(Ztm1X(:,end-2+1),[h,w]));
%of13estmator = opticalFlowLK('NoiseThreshold',0.009);estimateFlow(of13estmator, reshape(Ztm1X(:,end-3+1),[h,w]));
fgsOF = zeros(h,w,qX);
bgsOF = zeros(h,w,qX);
for t = 1 : qX
%t = 1; % Testing fame i = 1      
    fprintf('Testing fame %d at a measurement rate %2.2f \n', t, m/nX);
    fprintf('Updating proir by optical flow...\n');
    % Compute optical flow
    of12 = computeOF(reshape(Ztm1X(:,end-2+1),[h,w]), reshape(Ztm1X(:,end-1+1),[h,w]));
    of13 = computeOF(reshape(Ztm1X(:,end-3+1),[h,w]), reshape(Ztm1X(:,end-1+1),[h,w]));
    %of12 = estimateFlow(of12estmator, reshape(Ztm1X(:,end),[h,w]));
    %of13 = estimateFlow(of13estmator, reshape(Ztm1X(:,end),[h,w]));
    % Motion Compensation
    Ztm1X(:,end-2+1) = linearOFCompensate(reshape(Ztm1X(:,end-2+1),[h,w]), of12, 1,  true);
    Ztm1X(:,end-3+1) = linearOFCompensate(reshape(Ztm1X(:,end-3+1),[h,w]), of13, 0.5,true);
    
    % of = of(Ztm1X)
    % Performing CORPCA for separation
    fprintf('Performing CORPCA for separation...\n');
    ytX = PhiX*MX (:,t); % Input observation
    [xtX, vtX, ZtX, BtX] = corpca(ytX, PhiX, Ztm1X, Btm1X); % Performing CORPCA for separation
    % update prior information
    Ztm1X = ZtX;
    Btm1X = BtX;
    
    fgsOF(:,:,t) = reshape(xtX,[h,w]);
    bgsOF(:,:,t) = reshape(vtX,[h,w]);

    %fprintf('Recovered foreground error: %4.8f \n', norm(xtX - S(:,t),2)/(norm(S(:,t),2)));
    %fprintf('Recovered background error: %4.8f \n\n', norm(vtX - L(:,t),2)/(norm(L(:,t),2)));
    %reset(of12estmator);estimateFlow(of12estmator, reshape(Ztm1X(:,end-2+1),[h,w]));
    %reset(of13estmator);estimateFlow(of13estmator, reshape(Ztm1X(:,end-3+1),[h,w]));
end
toc

% ZtX(:,end) == vtX, which is the latest frame and prior of next frame

%%
%foreSave = VideoWriter('./fore_xing.avi');
%backSave = VideoWriter('./back_xing.avi');
%fore = reshape(xtX,[h,w]);
%imshow(reshape(vtX,[h,w]));
%imshow(reshape(Ztm1X(:,2),[h,w]));
%figure;montage(vRawClipGray(:,:,1:qX))
%figure;montage(uint8(abs(bgs)>0.1).*vRawClipGray(:,:,1:qX))
%figure;montage(uint8(abs(bgsOF)>0.1).*vRawClipGray(:,:,1:qX))
figure;montage(fgs.*reshape(MX,[h,w,qX]))
figure;montage(fgsOF.*reshape(MX,[h,w,qX]))

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