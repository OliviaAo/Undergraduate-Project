clear all
close all

bb=16;      % block size
RR=4;      % redundancy factor
K=RR*bb^2; % number of atoms in the dictionary
sigma=5; 

% pathForImages ='./test';
% imageName = '/1st_manual/01_manual1.gif';
pathForImages ='';
imageName = 'test3.png';
[IMin0,pp]=imread(strcat([pathForImages,imageName]));
IMin0=im2double(IMin0);
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;
end

IMin0 = mapminmax(IMin0, 0, 1); %对图片进行归一化处理
figure; imshow(IMin0,[]);

% FrangiScaleRange = [1 3];
% FrangiScaleRatio = 0.5;
% option.FrangiScaleRange = FrangiScaleRange;
% option.FrangiScaleRatio = FrangiScaleRatio;
% option.FrangiBetaOne = 0.5;
% option.FrangiBetaTwo = 8;
% option.BlackWhite = false;
% [IMin0,whatScale] = FrangiFilter2D(IMin0, option);
% figure; imshow(IMin0,[]);

% SegimageName = '/images/01_test.tif';
SegimageName = 'test3.3.png';
[IMin,pp]=imread(strcat([pathForImages,SegimageName]));
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

figure; imshow(IMin,[]);

FrangiScaleRange = [0.01 2];
FrangiScaleRatio = 0.05;
option.FrangiScaleRange = FrangiScaleRange;
option.FrangiScaleRatio = FrangiScaleRatio;
option.FrangiBetaOne = 0.6;
option.FrangiBetaTwo = 3;
option.BlackWhite = true;
[IMin,whatScale] = FrangiFilter2D(IMin, option);
figure; imshow(IMin,[]);
IMin=im2bw(IMin);
figure; imshow(IMin,[]);