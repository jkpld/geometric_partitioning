function [BW,S,I,iS] = createObjectImages(pixelList,Slist,Ilist,iSlist,numImRows)
% CREATEOBJECTIMAGES  Create the mask, edge, and intensity image for the
% object with pixels given by pixelList from an image with numImRows number
% of rows.

% James Kapaldo

% Get linear indices of object
j = rem(pixelList-1, numImRows) + 1; 
k = (pixelList - j)/numImRows + 1; 

r = [j,k];

% remove top left corner
minR = min(r)-1;
r = r-minR;

% Get the boundary size
maxR = max(r,[],1,'omitnan');

% Initialize mask
BW = false(maxR);
S = zeros(maxR);
I = zeros(maxR);
iS = zeros(maxR);

% Get linear indices for small image
inds = r(:,1) + (r(:,2)-1)*maxR(1);

% Assign in the values to the images.
BW(inds) = true;

if ~isempty(Slist)
    S(inds) = Slist;
end

if ~isempty(Ilist)
    I(inds) = Ilist;
end

if ~isempty(iSlist)
    iS(inds) = iSlist;
end

end
% createObjectImages changes log
% 2016-10-17 
% 2016-10-18/28 : undocumented changes
% 2016-10-29 : changed name to createObjectImages, renamed some inputs.