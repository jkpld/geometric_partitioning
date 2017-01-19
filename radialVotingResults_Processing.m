radialVotingResults_LD2P24 = radialVotingResults(:,1:3);
radialVotingResults_LD3P24 = radialVotingResults(:,(1:3)+1*4);
radialVotingResults_LD4P24 = radialVotingResults(:,(1:3)+2*4);
radialVotingResults_LD5P24 = radialVotingResults(:,(1:3)+3*4);
radialVotingResults_LD67P24 = radialVotingResults(:,(1:3)+4*4);


radialVotingResults_LD2P24(all(radialVotingResults_LD2P24==0,2),:) = [];
radialVotingResults_LD3P24(all(radialVotingResults_LD3P24==0,2),:) = [];
radialVotingResults_LD4P24(all(radialVotingResults_LD4P24==0,2),:) = [];
radialVotingResults_LD5P24(all(radialVotingResults_LD5P24==0,2),:) = [];
radialVotingResults_LD67P24(all(radialVotingResults_LD67P24==0,2),:) = [];


I = imread('testImage_image_LD4P24.tif');
objSize = round(size(I)/22);

idx = find(radialVotingResults_LD4P24(:,1)==0);
idx = [idx(:);size(radialVotingResults_LD4P24,1)];
offsets = 11;
for i = 2:numel(idx)-1
    radialVotingResults_LD4P24(idx(i):idx(i+1),2) = radialVotingResults_LD4P24(idx(i):idx(i+1),2) + offsets(i-1)*objSize(2);
end

I = imread('testImage_image_LD5P24.tif');
objSize = round(size(I)/22);

idx = find(radialVotingResults_LD5P24(:,1)==0);
idx = [idx(:);size(radialVotingResults_LD5P24,1)];
offsets = 11;
for i = 2:numel(idx)-1
    radialVotingResults_LD5P24(idx(i):idx(i+1),2) = radialVotingResults_LD5P24(idx(i):idx(i+1),2) + offsets(i-1)*objSize(2);
end


I = imread('testImage_image_LD67P24.tif');
objSize = round(size(I)/22);

idx = find(radialVotingResults_LD67P24(:,1)==0);
idx = [idx(:);size(radialVotingResults_LD67P24,1)];
offsets = [5,10,16];
for i = 2:numel(idx)-1
    radialVotingResults_LD67P24(idx(i):idx(i+1),2) = radialVotingResults_LD67P24(idx(i):idx(i+1),2) + offsets(i-1)*objSize(2);
end

radialVotingResults_LD2P24(:,1) = [];
radialVotingResults_LD3P24(:,1) = [];
radialVotingResults_LD4P24(:,1) = [];
radialVotingResults_LD5P24(:,1) = [];
radialVotingResults_LD67P24(:,1) = [];

rVR = {radialVotingResults_LD2P24, radialVotingResults_LD3P24, radialVotingResults_LD4P24, radialVotingResults_LD5P24, radialVotingResults_LD67P24};

%% Remove centers not from object

names = {'2','3','4','5','67'};
for i = 1:5
    BW = imread(['testImage_mask_LD' names{i} 'P24.tif']);
    I = imread(['testImage_image_LD' names{i} 'P24.tif']);
    mask = ~bwareaopen(I==0,20);
    BW = BW>0;
    BW = imdilate(BW,strel('disk',7));
    BW = double(BW);
    L = bwlabel(mask,8);
    
    % Combine points that are within 3 pixels of each other
    inds = getPdistInds(size(rVR{i},1));
    d = pdist(rVR{i});
    toCombine = find(d<3);
    while ~isempty(toCombine)
    
        inds = inds(toCombine,:);
        rVR{i}(inds(:,1),:) = (rVR{i}(inds(:,1),:) + rVR{i}(inds(:,2),:))/2;
        rVR{i}(inds(:,2),:) = [];
        
        inds = getPdistInds(size(rVR{i},1));
        d = pdist(rVR{i});
        toCombine = find(d<3);
    end
    
    in = interp2(BW,rVR{i}(:,1),rVR{i}(:,2));
    rVR{i}(in<1,:) = [];
    objNum = interp2(L,rVR{i}(:,1),rVR{i}(:,2));
    rVR{i} = [objNum, rVR{i}];
    rVR{i}(objNum==0 | isnan(objNum),:) = [];    
end


%%
BW = imread('testImage_mask_LD67P24.tif');
BW = BW>0;

figure
imshow(BW)
hold on
plot(rVR{5}(:,2),rVR{5}(:,3),'.r')

%%

NS = createns(markedCenters_LD3P24(:,[3,2]));

tc = rVR{2}(:,2:3);
objNum = rVR{2}(:,1);
idx5 = rangesearch(NS,tc,5);
idx3 = rangesearch(NS,tc,3);

F5_rVR = sum(~cellfun(@isempty,idx5))/size(tc,1);
F3_rVR = sum(~cellfun(@isempty,idx3))/size(tc,1);

% Get the number difference between the number of centers found and
% the correct number of centers
disp here
numCents = accumarray(objNum,1,[],[],0);
disp here
correctNumCents = accumarray(markedCenters_LD3P24(:,1),1,[],[],0);
objNumbers = unique(markedCenters_LD3P24(:,1));

d = numCents(objNumbers) - correctNumCents(objNumbers);
d(abs(d)>3) = [];
d = d + 4;
N_rVR = permute(accumarray(d,1,[7,1],[],0),[3,2,1]);


