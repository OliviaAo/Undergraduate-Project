clear all
close all

pathForImages ='./test';
imageName = '/images/01_test.tif';
[IMin0,pp]=imread(strcat([pathForImages,imageName]));
IMin0=im2double(IMin0);
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;
end

SegimageName = '/1st_manual/01_manual1.gif';
[IMin,pp]=imread(strcat([pathForImages,SegimageName]));
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

index = find ( IMin );
vessel = zeros( size(IMin0));
vessel(index) = IMin0(index);
figure; imshow( vessel,[]);

FrangiScaleRange = [0.1 3];
FrangiScaleRatio = 0.5;
option.FrangiScaleRange = FrangiScaleRange;
option.FrangiScaleRatio = FrangiScaleRatio;
option.FrangiBetaOne = 0.5;
option.FrangiBetaTwo = 8;
% option.BlackWhite = false;
[vessel,whatScale] = FrangiFilter2D(vessel, option);
figure; imshow(vessel,[]);


H = fspecial('disk',2);
blurred = imfilter(IMin0,H,'replicate');
index = find (~IMin);
vessel(index) = blurred(index);

figure; imshow( vessel,[]);

FrangiScaleRange = [0.1 3];
FrangiScaleRatio = 0.5;
option.FrangiScaleRange = FrangiScaleRange;
option.FrangiScaleRatio = FrangiScaleRatio;
option.FrangiBetaOne = 0.5;
option.FrangiBetaTwo = 8;
% option.BlackWhite = false;
[vessel,whatScale] = FrangiFilter2D(vessel, option);
figure; imshow(vessel,[]);
