bw = imread('K:\Google_Drive\MATLAB\seed_point_based_segmentation\exampleImages\testNuclei_mask.tif');

figure(1)
imshow(bw)

imSize = size(bw);
numObjs = 1;

options = partitionOptions;
objectScale = double(max(reshape(bwdist(~bw),[],1)));

[B,n,kappa] = computeBoundaryInformation(bw,objectScale,options);

locs = max_curvature_idx(kappa);


figure(2)
imagesc(tmpI)
daspect([1 1 1])
colorbar

tmpI = zeros(size(bw)); %#ok<*UNRCH>
for i = 1:numObjs
    idx = B{i}(1:end-1,1) + (B{i}(1:end-1,2)-1)*imSize(1);
    good = ~isnan(idx);
    tmpI(idx(good)) = kappa{i}(good);
    
    line(B{i}(locs{i},2),B{i}(locs{i},1),'marker','x','color','g','linestyle','none');
    
end



cmap = createColorGradient([0 0 1],[0 0 0],10);
cmap = [cmap; createColorGradient([0 0 0], [1 0 0],10)];
colormap(cmap)
clim([-0.1,0.1])
drawnow;
setTheme(gcf,'dark');


%%

T = [n{1}(:,2),-n{1}(:,1)];
edgeInds = delaunayDeclumping(B{1},T,locs{1},imSize,1,1);

bw2 = bw;
bw2(edgeInds) = 0;

figure(4)
imshow(bw2)

%%

pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];

ind = 5;
names = {'2','3','4','5','67'};
I = imread(im_pth(names{ind}));
BW = imread(bw_pth(names{ind})) > 0;
L = bwlabel(imfill(I~=0,'holes'));

results = load(results_pth(names{ind}));
result_fields = fieldnames(results);
tmp = results.(result_fields{1});
truthDat = tmp;
r0x = accumarray(truthDat(:,1), truthDat(:,2), [484,1],@(x) {x});
r0y = accumarray(truthDat(:,1), truthDat(:,3), [484,1],@(x) {x});
r0 = cellfun(@(x,y) [x,y], r0x, r0y, 'UniformOutput',false);

options = partitionOptions();
options.Use_GPU = true;
[BW_part, cuts, Info] = partition_objects(I,BW,r0,options);

figure
imshow(BW_part)
hold on
plot(truthDat(:,3),truthDat(:,2),'.r')