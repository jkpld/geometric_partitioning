[I,~,BW] = getDeclumpTestCase('6and7','Cancer');
CC = bwconncomp(BW,8);

%%

labels = labelNucleiCenters2(I,BW,CC);





%%

imsize = round(size(I)/22);

[j,k] = ind2sub([22 22],205);
% sI = I((j-1)*imsize(1) + (1:97)+1, (k-1)*imsize(2) + (1:92)+1);
% sBW = BW((j-1)*imsize(1) + (1:97)+1, (k-1)*imsize(2) + (1:92)+1);
% sI = I((j-1)*imsize(1) + (1:114)+1, (k-1)*imsize(2) + (1:99)+1);
% sBW = BW((j-1)*imsize(1) + (1:114)+1, (k-1)*imsize(2) + (1:99)+1);
sI = I((j-1)*imsize(1) + (1:114)+1, (k-1)*imsize(2) + (1:94)+1);
sBW = BW((j-1)*imsize(1) + (1:114)+1, (k-1)*imsize(2) + (1:94)+1);

B = bwboundaries(sBW);
B = B{1};

sI = double(sI);

% PAD_SIZE = 10;
% sIp = padarray(sI,[PAD_SIZE,PAD_SIZE],'symmetric');
% sBWp = padarray(sBW,[PAD_SIZE,PAD_SIZE],0);

% sBW = ~bwareaopen(~sBW,options.Minimum_Hole_Size,4);
% CC = bwconncomp(sBW,8);

% figure
% imshow(sI,[])
options = declumpOptions();
options.Potential_Depth = -1;
options.Potential_Minimum_Location = 2;
options.Potential_Extent = 15;
options.Use_GPU = true;
options.Debug = true;
[declumpedBW, cuts, Info] = declumpNuclei(sI,sBW,options);


sI = imfilter(sI,fspecial('gaussian',7,1));

%%

H = fspecial('disk',5);


tmpI = sIp./imfilter(sIp,H,'symmetric');
minMaxI = prctile(tmpI(:),[10,90]);
tmpI = tmpI - minMaxI(1);
tmpI = tmpI/(minMaxI(2)-minMaxI(1));
tmpI(tmpI<0) = 0;
tmpI(tmpI>1) = 1;
% tmpI(~imerode(sBWp,strel('disk',3))) = 1;
figure(10)
clf(10)
imshow(tmpI,[])
colorbar
clim([0 1])

Z = @(X,Y,r) cos(pi/r * X).*cos(pi/r * Y).*exp(-(X).^2/(2*7^2)-(Y).^2/(2*7^2))/numel(X);

[X,Y] = meshgrid(-15:15);
[Th,R] = cart2pol(X,Y);
rXY = [cosd(45), -sind(45); sind(45), cosd(45)]*[X(:),Y(:)]';
rX = reshape(rXY(1,:),size(X));
rY = reshape(rXY(2,:),size(Y));

ap = ones(size(sIp));
abw = ones(size(sIp));
% meanI = mean(tmpI(logical(sBWp)));
% tmpI = tmpI - meanI;
for i = 5:10
H = Z(X,Y,i);
rH = Z(rX,rY,i);

tap = abs(imfilter(sIp,H,'symmetric')) + abs(imfilter(sIp,rH,'symmetric'));
% tap = tap.*abs(imfilter(sIp,H,'symmetric')) + abs(imfilter(sIp,rH,'symmetric'));
% tbw = abs(imfilter(double(sBWp),H,'symmetric')) + abs(imfilter(double(sBWp),rH,'symmetric'));
% ap = ap.*imfilter(tap,fspecial('gaussian',35,5));%medfilt2(tap,11*[1 1],'symmetric');%
% abw = abw.*imfilter(tbw,fspecial('gaussian',35,5));%medfilt2(tbw,11*[1 1],'symmetric');%
% ap = ap.*abs(imfilter(tap,fspecial('gaussian',35,5)) - imfilter(tbw,fspecial('gaussian',35,5)));
ap = ap.*tap;

end
% tmpI = tmpI + meanI;
ap = imfilter(ap,fspecial('gaussian',37,5),'symmetric');
ap = ap/max(ap(:));
% ap = ap./(abw+1);
% ap = medfilt2(ap,15*[1 1],'symmetric');
% ap = imdilate(ap,strel('disk',5));
ap = ap/max(ap(:));


% sIpf = imfilter(tmpI,fspecial('gaussian',35,5),'symmetric');
sIpf = imfilter(tmpI,fspecial('disk',5),'symmetric');
tmpI = tmpI.*(1-ap) + ap;
tmpI(tmpI<0.05) = 0.05;
tmpI = imfilter(tmpI,fspecial('gaussian',7,1),'symmetric');
% tmpI = imopen(tmpI,strel('disk',4));
% tmpI = imfilter(tmpI,fspecial('disk',3),'symmetric');
% tmpI = tmpI*0.8 + 0.2;

% H = fspecial('disk',11);
% tmpI = tmpI./imfilter(tmpI,H,'symmetric');
% tmpI = tmpI-min(tmpI(:)) + 0.05;
% tmpI = tmpI/max(tmpI(:));

abw = abw(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
ap = ap(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
tmpI = tmpI(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);

figure(16)
clf(16)
imshow(tmpI,[])
line(B(:,2),B(:,1),'color','r')
colorbar
clim([0 1])



figure(11)
clf(11)
imshow(ap,[])
line(B(:,2),B(:,1),'color','r')
% colorbar

figure(12)
clf(12)
imagesc(Info{1}.V.*(1./tmpI))
line(B(:,2),B(:,1),'color','r')
daspect([1 1 1])
colorbar
clim([0,1])

figure(17)
clf(17)
imagesc(abw)
daspect([1 1 1])
colorbar
% clim([0,1])

%%


% currentmean = startmean;
% for i = 1 : ng
%     prevmean = currentmean;
%     ero = imerode(ero,strel('diamond',1));
%     rec = imreconstruct(ero,img,4);
%     currentmean = mean2(rec);
%     gs(i) = prevmean - currentmean;
% end
% gs = gs .* (100 / startmean);

gS = 1:15;
gI = zeros([size(sI),numel(gS)+1]);
gI(:,:,1) = sI;
ero = imerode(sI,strel('disk',11));
se = strel('disk',1);
gI(:,:,2) = imreconstruct(ero,sI,4);
% for i = 2:numel(gS)
%     ero = imerode(ero,se);
%     if i == numel(gS)
%         figure
%         imshow(ero,[])
%         colorbar
%     end
%     gI(:,:,i+1) = imreconstruct(ero,sI,4);
% end
% gI = imfilter(gI,fspecial('disk',4));
% gI = flip(diff(flip(gI,3),1,3),3);
figure
imshow(gI(:,:,2),[])


%%
% I = double(I);
tmpI = single(I);
time0 = tic;

tmpI_d = gpuArray(tmpI);
tmpI_d = imfilter(tmpI_d,fspecial('gaussian',7,1),'symmetric');
tmpI = gather(tmpI_d);
% clear tmpI_d
tmpI_d = imreconstruct(gpuArray(imerode(tmpI,strel('disk',7))),tmpI_d,4);
tmpI = gather(tmpI_d);
clear tmpI_d
tmpI = imopen(tmpI,strel('disk',5));
tmpI_d = gpuArray(tmpI);
% figure(20)
% clf(20)
% imshow(tmpI,[])
% H = fspecial('disk',4);
H = fspecial('gaussian',21,3);
tmpI_d = tmpI_d./imfilter(tmpI_d,H,'symmetric');
% tmpI = imfill(tmpI);
tmpI = gather(tmpI_d);
clear tmpI_d
CC = bwconncomp(BW,8);
minI = cellfun(@(x) prctile(tmpI(x),10),CC.PixelIdxList);
maxI = cellfun(@(x) prctile(tmpI(x),30),CC.PixelIdxList);
% minI = double(minI);
% maxI = double(maxI);
time1 = tic;
for i = 1:max(L(:))
%     mask = L==i;
%     minMaxI = prctile(tmpI(mask),[2,30]);
    tmpI(CC.PixelIdxList{i}) = tmpI(CC.PixelIdxList{i}) - minI(i);
    tmpI(CC.PixelIdxList{i}) = tmpI(CC.PixelIdxList{i})/(maxI(i)-minI(i));
end
toc(time1)
% minMaxI = prctile(tmpI(:),[10,60]);
% tmpI = tmpI - minMaxI(1);
% tmpI = tmpI/(minMaxI(2)-minMaxI(1));
tmpI(tmpI<0) = 0;
tmpI(tmpI>1) = 1;

tmpI = tmpI*0.999 + 0.001;
tmpI = imfilter(tmpI,fspecial('gaussian',7,1),'symmetric');


tmpI = tmpI.*single(BW);
tmpI(~isfinite(tmpI)) = 0;
tmpI = imfill(tmpI);
tmpI = tmpI.*single(BW);
tmpI = gather(tmpI);


toc(time0)

% figure(20)
% clf(20)
% imshow(sI,[])

%%
figure(21)
clf(21)

imshow(1./tmpI,[])
%%
figure(22)
clf(22)
imagesc(Info{1}.V.*(1./tmpI))
line(B(:,2),B(:,1),'color','r')
daspect([1 1 1])
colorbar
clim([0,1])

%%

figure(23)
clf(23)
imagesc(Info{1}.V)
line(B(:,2),B(:,1),'color','r')
daspect([1 1 1])
colorbar
clim([0,1])

%%

tmpI = imreconstruct(imerode(sI,strel('disk',11)),sI,4);
tmpI = imopen(tmpI,strel('disk',5));
tmpI = sI./tmpI;%imreconstruct(imerode(sI,strel('disk',7)),sI,4);
tmpI(~imerode(sBW,strel('disk',3))) = 1;

mu = imfilter(tmpI,fspecial('disk',7),'symmetric');
Istd = sqrt(imfilter((tmpI-mu).^2,fspecial('disk',7),'symmetric'));
Istd_mask = imerode(Istd>0.2,strel('disk',3)) & sBW;

% Istd = imfilter(Istd,fspecial('gaussian',35,5),'symmetric');
% Istd = medfilt2(Istd,4*[1 1],'symmetric');
figure
imshow(Istd,[])
line(B(:,2),B(:,1),'color','r')
B2 = bwboundaries(Istd_mask);
for i = 1:numel(B2)
line(B2{i}(:,2),B2{i}(:,1),'color','b')
end