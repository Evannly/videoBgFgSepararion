function imgComp = linearOFCompensate(img, ofObj, factor, vectorize)
% Make linear motion compensation to img according to ofObj

% Inputs:
%   img: original image, should be gray scale image
%   ofObj: optical flow struct obtain by Matlab offical CV toolbox, should
%       have the format:
%             flow = 
%   
%               opticalFlow with properties:
%   
%                          Vx: [60×80 single]
%                          Vy: [60×80 single]
%                 Orientation: [60×80 single]
%                   Magnitude: [60×80 single]
%   factor: scale factor of compensation, i.e. x' = factor*motion + x.
%       Default to be 1
%   vecize: if vectorize result or not

% Outputs:
%   imgComp: compensated image.
    if nargin == 2
        factor = 1;
        vectorize = true;
    elseif nargin == 3
        vectorize = true;
    end
    
    imgComp = zeros(size(img));
    [h, w, c] = size(img);
    % x'(yy,xx) = x(y,x), 
    % where yy = y + factor*motiony, xx = factor+motionx
    for y = 1:h
        for x = 1:w
            %yy = y + factor * bound(round(ofObj.Vy(y,x)), 1, h);
            %xx = x + factor * bound(round(ofObj.Vx(y,x)), 1, w);
            yy = y + round(factor * ofObj.Vy(y,x));
            xx = x + round(factor * ofObj.Vx(y,x));
            if (yy>=1) && (yy<=h) && (xx>=1) && (xx<=w)
                imgComp(yy, xx, :) = img(y, x, :);
            end
        end
    end
    if vectorize
        imgComp = reshape(imgComp,[h*w*c, 1]);
    end
end

function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end