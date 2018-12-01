classdef utils

   
   methods(Static)
        %% util functions
        function tansferedImgs = getTansferedImgs(fpath,c,nClip)
            resizeScale = 0.2;
            imagefiles = dir([fpath '*.ppm']);
            if nargin == 2
                nfiles = length(imagefiles); % Number of files found
            else
                nfiles = nClip;
            end
            
            c0 = 1;
            [h0, w0, c0] = size(imread([fpath, imagefiles(1).name]));

            h = floor(h0 /resizeScale);
            w = floor(w0 /resizeScale);
            c = min(c, c0);          % Correct chanUse
            %images = zeros([h*w,nfiles]);                                  % image size (120,160,3)
            tansferedImgs = zeros([h , w , c, nfiles]); % image size (120,160,3)

            for ii = 1:nfiles
                currentfilename = imagefiles(ii).name;
                currentimage = imread([fpath, currentfilename]);

                if c == 1
                    currentimage = imresize(rgb2gray(currentimage), 1/resizeScale);
                else
                    currentimage = imresize(currentimage, 1/resizeScale);
                end

                tansferedImgs(:,:,:, ii) = currentimage;
            end

        end 

        
                %% Utility functions
        function of = computeOF(img1, img2, method)

            if nargin == 2
                method = 'LK';
            end

            if strcmp(method, 'LK')
                ofestmator = opticalFlowLK('NoiseThreshold', 0.009);
                estimateFlow(ofestmator, img1);
                of = estimateFlow(ofestmator, img2);
            end

        end

        function saveVideo(vFileName, vVar, frameRate, vFileSuffix)
            if nargin ==1
                error('Too few inputs!');
            elseif nargin == 2
                frameRate = 30;
                vFileSuffix = 'avi';
                vFileFormat = 'Motion JPEG AVI';
            elseif nargin ==3
                vFileSuffix = 'avi';
                vFileFormat = 'Motion JPEG AVI';
            else
                if strcmp(vFileSuffix, 'avi')
                    vFileFormat = 'Motion JPEG AVI';
                elseif strcmp(vFileSuffix, 'mp4')
                    vFileFormat = 'MPEG-4';
                end
            end



            vFileName = [vFileName '.' vFileSuffix];

            if isa(vVar, 'uint8')
                vVar = double(vVar);        
            end

            vVar = (vVar - min(vVar(:))) / max(vVar(:));
            disp(vFileFormat)
            vFile = VideoWriter(vFileName, vFileFormat);
            vFile.FrameRate = frameRate;
            open(vFile);
            writeVideo(vFile, vVar);
            close(vFile);
        end

   end
end