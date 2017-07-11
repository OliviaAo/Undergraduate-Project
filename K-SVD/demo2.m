%============================================================
%               demo2 - denoise an image
% this is a run_file the demonstrate how to denoise an image, 
% using dictionaries. The methods implemented here are the same
% one as described in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
%============================================================
clear all
close all

bb=8;      % block size
RR=8;      % redundancy factor
K=RR*bb^2; % number of atoms in the dictionary
sigma=5; 
% Rect = [286  236  31  31];
Rect_std = [136  126  31  31];
% Rect_std = [136  126  7  7];
% Rect_std1 = [185 435 7 7]
pathForImages ='./test';
imageName = '/1st_manual/01_manual1.gif';
segI = mapminmax( im2double(imread(strcat([pathForImages,imageName]))),0, 1);

[NN1,NN2] = size(segI);
IMout=zeros(NN1,NN2);
BB=32;
col = 136; row = 126;
IMout(row:row+BB-1,col:col+BB-1)=segI(row:row+BB-1,col:col+BB-1);
figure;
subplot(1,3,1); imshow(segI,[]);
subplot(1,3,2); imshow(IMout,[]);

[IMin0,pp] = imcrop ( imread(strcat([pathForImages,imageName])), Rect_std);
IMin0=im2double(IMin0);
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;
end
IMin0 = mapminmax(IMin0, 0, 1); %????????????????????
figure; 
imshow(IMin0,[]);

SegimageName = '/images/01_test.tif';
[IMin] =rgb2gray(imread(strcat([pathForImages,SegimageName])));

IMin  = imadjust( IMin );
figure;
subplot(1,3,1); imshow(IMin,[]);  title('before histeq');

IMin  = adapthisteq( IMin, 'NumTiles',[4 4],'ClipLimit',0.01);
subplot(1,3,2); imshow(IMin,[]); title('after histeq');

IMin= double(IMin);

IMin_back = IMin;
FrangiScaleRange = [0.1 5];
FrangiScaleRatio = 0.2;
option.FrangiScaleRange = FrangiScaleRange;
option.FrangiScaleRatio = FrangiScaleRatio;
option.FrangiBetaOne = 0.4;
option.FrangiBetaTwo = 25;
option.BlackWhite = true;
[IMin_back,whatScale] = FrangiFilter2D(IMin_back, option);


% Get blocks from IMin:
BB = 128;
slidingDis = 32;
[blocks_Rect,idx1] = my_im2col(IMin,[BB,BB],slidingDis);
[NN1,NN2] = size(IMin);
[N1,N2]=size(blocks_Rect);
[rows,cols] = ind2sub(size(IMin)-BB+1,idx1);
for i=1:N2
    block =reshape(blocks_Rect(:,i),[BB,BB]);
    s = std(blocks_Rect(:,i),0,1);
    disp(['??: ',num2str(s)]);
    col = cols(i); row = rows(i);
    IMout = zeros(NN1,NN2);
    IMout(row:row+BB-1,col:col+BB-1)=IMout(row:row+BB-1,col:col+BB-1)+block;
%     figure;imshow(IMout,[]);
    if(s>50)
        FrangiScaleRange = [0.1 5];
        FrangiScaleRatio = 0.2;
        option.FrangiScaleRange = FrangiScaleRange;
        option.FrangiScaleRatio = FrangiScaleRatio;
        option.FrangiBetaOne = 0.4;
        option.FrangiBetaTwo = 40;
        option.BlackWhite = true;
        [block,whatScale] = FrangiFilter2D(block, option);
        blocks_Rect(:,i) = im2col(block,[BB BB],'distinct');
%         figure;imshow(block,[]);
    else
        FrangiScaleRange = [0.1 5];
        FrangiScaleRatio = 0.2;
        option.FrangiScaleRange = FrangiScaleRange;
        option.FrangiScaleRatio = FrangiScaleRatio;
        option.FrangiBetaOne = 0.4;
        option.FrangiBetaTwo = 45;
        option.BlackWhite = true;
        [block,whatScale] = FrangiFilter2D(block, option);
        blocks_Rect(:,i) = im2col(block,[BB BB],'distinct');
%         figure;imshow(block,[]);
%         block =reshape(blocks_Rect(:,i),[BB,BB]);
%         figure;imshow(block,[]);
    end    
end

[NN1,NN2] = size(IMin);
Weight=zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
for i=1:N2
%     disp(['ok1']);
    block =reshape(blocks_Rect(:,i),[BB,BB]);
%     figure; imshow(block,[]);
%     [IoutGlobal,output] = denoiseImageGlobal(D, block,sigma,K);
%     figure; imshow(IoutGlobal,[]);
    col = cols(i); row = rows(i);
    IMout(row:row+BB-1,col:col+BB-1)=IMout(row:row+BB-1,col:col+BB-1)+block;
    Weight(row:row+BB-1,col:col+BB-1)=Weight(row:row+BB-1,col:col+BB-1)+ones(BB);
%     disp(['ok2']);
end
index  = find (Weight);
IOut  = zeros ( size (IMout));
IOut(index) = IMout(index)./Weight(index);
subplot(1,3,3); imshow(IOut,[]); title('after frangi');

IMin = IOut;


% Get dictionary based on IMin0
slidingDis = 5;
[blocks,idx] = my_im2col(segI,[bb,bb],slidingDis);
D = selectDic ( blocks, K );
vecOfDic = mean(D);
D = D - repmat(vecOfDic,size(D,1),1);

figure;
I = displayDictionaryElementsAsImage(D, floor(sqrt(K)), floor(size(D,2)/floor(sqrt(K))),bb,bb);
title('The dictionary trained on patches from segmented images');

%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   O V E R C O M P L E T E 
%                        D C T    D I C T I O N A R Y
% ==========================================================================

[NN1,NN2] = size(IMin);
Weight=zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(IMin)-BB+1,idx1);
for i=1:N2
    block =reshape(blocks_Rect(:,i),[BB,BB]);
    [IoutDCT,output] = denoiseImageDCT(block,sigma,K);
    col = cols(i); row = rows(i);
    IMout(row:row+BB-1,col:col+BB-1)=IMout(row:row+BB-1,col:col+BB-1)+IoutDCT;
    Weight(row:row+BB-1,col:col+BB-1)=Weight(row:row+BB-1,col:col+BB-1)+ones(BB);
end

IOut  = zeros ( size (IMout));
IOut = (IMin+0.034*sigma*IMout)./(1+0.034*sigma*Weight);

% pathForImages ='./test';
% image1 = '/1st_manual/10_manual1.gif';
% [Image1,pp] =imread(strcat([pathForImages,image1]));
Image1 = IMin;

[m,n] = size(IOut);
A = double(Image1);
B = double(IOut);
C = sum( sum( (A-B).^2 ) );
MSE = C / (m * n);
PSNR = 10*log10( (255^2) / MSE );                                                        
disp(['DCT psnr:',num2str(PSNR)]);

figure;
subplot(1,3,1); imshow(IOut,[]); title('Clean Image by DCT Trained dictionary');
% ==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   G L O B A L 
%           ( O R   G I V E N )   D I C T I O N A R Y
%==========================================================================

[NN1,NN2] = size(IMin);
Weight=zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(IMin)-BB+1,idx1);
for i=1:N2
    block =reshape(blocks_Rect(:,i),[BB,BB]);
    [IoutGlobal,output] = denoiseImageGlobal(D, block,sigma,K);
    col = cols(i); row = rows(i);
    IMout(row:row+BB-1,col:col+BB-1)=IMout(row:row+BB-1,col:col+BB-1)+IoutGlobal;
    Weight(row:row+BB-1,col:col+BB-1)=Weight(row:row+BB-1,col:col+BB-1)+ones(BB);
end

IOut1  = zeros ( size (IMout));
IOut1 = (IMin+0.034*sigma*IMout)./(1+0.034*sigma*Weight);

[m,n] = size(IOut1);
A = double(Image1);
B = double(IOut1);
C = sum( sum( (A-B).^2 ) );
MSE = C / (m * n);
PSNR = 10*log10( (255^2) / MSE );                                                        
disp(['Global psnr:',num2str(PSNR)]);
subplot(1,3,2); imshow(IOut1,[]); title('Clean Image by Global Trained dictionary');

%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   A   D I C T  I O N A R Y
%                  T R A I N E D   O N   N O I S Y   I M A G E
%==========================================================================

[NN1,NN2] = size(IMin);
Weight=zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(IMin)-BB+1,idx1);
for i=1:N2
    disp(['ok1']);
    block =reshape(blocks_Rect(:,i),[BB,BB]);
%     figure; imshow(block,[]);
    [IoutAdaptive,output] = denoiseImageKSVD(D, block,sigma,K);
%     figure; imshow(IoutGlobal,[]);
    col = cols(i); row = rows(i);
    IMout(row:row+BB-1,col:col+BB-1)=IMout(row:row+BB-1,col:col+BB-1)+IoutAdaptive;
    Weight(row:row+BB-1,col:col+BB-1)=Weight(row:row+BB-1,col:col+BB-1)+ones(BB);
    disp(['ok2']);
end

IOut2  = zeros ( size (IMout));
IOut2 = (IMin+0.034*sigma*IMout)./(1+0.034*sigma*Weight);

[m,n] = size(IOut2);
A = double(Image1);
B = double(IOut2);
C = sum( sum( (A-B).^2 ) );
MSE = C / (m * n);
PSNR = 10*log10( (255^2) / MSE );                                                        
disp(['KSVD psnr:',num2str(PSNR)]);
subplot(1,3,3); imshow(IOut2,[]); title('Clean Image by Adaptive  dictionary');

figure;
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
title('The dictionary trained on patches from segmented images');
