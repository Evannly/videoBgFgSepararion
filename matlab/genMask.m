% %% generate the mask 
% % for the foreground
% % figure
% % subplot(2,2,1);hist(fgs(:));title('fgs');subplot(2,2,2);hist(bgs(:));title('bgs');
% % subplot(2,2,3);hist(fgsOF(:));title('fgsOF');subplot(2,2,4);hist(bgsOF(:));title('fgsOF');
% % message = sprintf('fgs: mean = %d , median = %d, max = %d, min = %d ', mean(fgs(:)),median(fgs(:)), max(fgs(:)), min(fgs(:)));
% % disp(message);
% 
% % encode with Uint8
% % fgs = uint8(fgs);bgs = uint8(bgs);
% % fgsOF = uint8(fgsOF);bgsOF = uint8(bgsOF);
% 
% 
% frameRate = 15;
% u = utils;
% %% corpca
% path = '/home/exxact/Documents/xiaojian/projects/DeepStyleTransfer/video_output/';
% channel = 3;
% % generate mask
% thershold =1 ;
% fgsMsk = fgs;
% fgsMsk(fgsMsk < thershold) = 0;
% fgsMsk = sign(fgsMsk);
% bgsMsk = 1 - fgsMsk;
% % sythnethe the video	
% bgsTans = u.getTansferedImgs(path,channel);
% sythVideo = fgsMsk.* fgs + bgsMsk .* bgsTans;
% u.saveVideo('./fbVideos/v1/sythVideo', abs(sythVideo),frameRate);
% 
% %% corpca-OF
% path = '/home/exxact/Documents/xiaojian/projects/DeepStyleTransfer/video_output1/'
% channel = 3;
% % generate mask
% thershold =0 ;
% fgsOFMsk = fgsOF;
% fgsOFMsk(fgsOFMsk < thershold) = 0;
% fgsOFMsk = sign(fgsOFMsk);
% bgsOFMsk = 1 - fgsOFMsk;
% % sythnethe the video	
% bgsTans = u.getTansferedImgs(path,channel);
% sythOFVideo = fgsOFMsk.* fgsOF + bgsOFMsk .* fgsOFTrans;
% u.saveVideo('./fbVideos/sythOFVideo', abs(sythOFVideo),frameRate);
% %%
% % readVideo to mat
% nClip = 100;c0 = 1;c = 3;resizeScale = 1;
% imagefiles = VideoReader('./videos/bootstrap.avi');
% vRaw = read(imagefiles, [1 nClip]); % Read video. Size: (60,80,3,nClip)
% u.saveVideo('./fbVideos/oriVideoClip', abs(vRaw),frameRate);
% 
% % Read iamges to mat
% path = './videos/bootstrap/';
% channel = 3; nClip = 100;
% oriVideoClip = u.getTansferedImgs(path,channel,nClip);
% u.saveVideo('./fbVideos/oriVideoClip', abs(oriVideoClip),frameRate);
% 
% deal with the self imges
path = './videos/bridge/';
channel = 1; nClip = 30;
oriVideoClip = u.getTansferedImgs(path,channel,nClip);
u.saveVideo('./fbVideos/bri/oriBridgeVideoClip', abs(oriVideoClip),frameRate);

path = '/home/exxact/Documents/xiaojian/projects/DeepStyleTransfer/video_output/';
channel = 3;resizeScale = 1;
% generate mask
thershold =5 ;
fgsMsk = fgs;
fgsMsk(fgsMsk < thershold) = 0;
fgsMsk = sign(fgsMsk);
bgsMsk = 1 - fgsMsk;
% sythnethe the video	
bgsTans = u.getTansferedImgs(path,channel);
sythVideo = uint8(fgsMsk.* fgs + bgsMsk .* bgsTans);
u.saveVideo(['./fbVideos/bri/sythVideo_' num2str(thershold)], abs(sythVideo),frameRate);


