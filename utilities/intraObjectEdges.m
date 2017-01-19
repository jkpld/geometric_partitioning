function S = intraObjectEdges(I,BW,CC,options)
% INTRAOBJECTEDGES  Get the edges between sub-objects in an object
%
% S = intraObjectEdges(I,BW,CC,options)
%
% I - input image
% BW - object mask
% CC - connected component structure for BW
% options - declumpOptions object
%
% S - image with the same size as I and class single. regions of inbetween
% sub-objects have a value > 1, all other values are 1.

% James Kapaldo
% 2016-11-28

params.EROSION_SIZE = 7;
params.LOWER_EDGE_THRESHOLD_PERCENT = 10;
params.UPPER_EDGE_THRESHOLD_PERCENT = 30;

if options.Use_GPU
    S = intraObjectEdges_GPU(I,BW,CC,params);
else
    S = intraObjectEdges_CPU(I,BW,CC,params);
end

end

function S = intraObjectEdges_CPU(I,BW,CC,params)

EROSION_SIZE = params.EROSION_SIZE;
THRESH = [params.LOWER_EDGE_THRESHOLD_PERCENT, params.UPPER_EDGE_THRESHOLD_PERCENT];

S = single(I);

S = imfilter(S,fspecial('gaussian',7,1),'symmetric');
S = imreconstruct(gpuArray(imerode(S,strel('disk',EROSION_SIZE))),S,4);
S = imopen(S,strel('disk',5));

H = fspecial('gaussian',21,3);
S = S./imfilter(S,H,'symmetric');

minI = cellfun(@(x) prctile(S(x),THRESH(1)),CC.PixelIdxList);
maxI = cellfun(@(x) prctile(S(x),THRESH(2)),CC.PixelIdxList);

for i = 1:CC.NumObjects
    S(CC.PixelIdxList{i}) = (S(CC.PixelIdxList{i}) - minI(i)) / (maxI(i)-minI(i));
end

S(S<0) = 0;
S(S>1) = 1;

S = S*0.999 + 0.001;
S = imfilter(S,fspecial('gaussian',7,1),'symmetric');

S = S.*single(BW);
S(~isfinite(S)) = 0;
S = imfill(S);
S = (1./S);
S(~BW) = 1;


end

function S = intraObjectEdges_GPU(I,BW,CC,params)

EROSION_SIZE = params.EROSION_SIZE;
THRESH = [params.LOWER_EDGE_THRESHOLD_PERCENT, params.UPPER_EDGE_THRESHOLD_PERCENT];

S = single(I);

S_d = gpuArray(S);
S_d = imfilter(S_d,fspecial('gaussian',7,1),'symmetric');
S = gather(S_d);

S_d = imreconstruct(gpuArray(imerode(S,strel('disk',EROSION_SIZE))),S_d,4);
S = gather(S_d);
clear S_d

S = imopen(S,strel('disk',5));
S_d = gpuArray(S);

H = fspecial('gaussian',21,3);
S_d = S_d./imfilter(S_d,H,'symmetric');

S = gather(S_d);
clear S_d

minI = cellfun(@(x) prctile(S(x),THRESH(1)),CC.PixelIdxList);
maxI = cellfun(@(x) prctile(S(x),THRESH(2)),CC.PixelIdxList);


for i = 1:CC.NumObjects
    S(CC.PixelIdxList{i}) = S(CC.PixelIdxList{i}) - minI(i);
    S(CC.PixelIdxList{i}) = S(CC.PixelIdxList{i})/(maxI(i)-minI(i));
end

S(S<0) = 0;
S(S>1) = 1;

S = S*0.999 + 0.001;
S = imfilter(S,fspecial('gaussian',7,1),'symmetric');

S = S.*single(BW);
S(~isfinite(S)) = 0;
S = imfill(S);
S = (1./S);
S(~BW) = 1;

end