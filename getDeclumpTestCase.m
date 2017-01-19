function [I,Is,BW] = getDeclumpTestCase(caseNum,cellType)

if nargin == 1
    switch caseNum
        case 1 % Nuclei clump with 8 nuclei, 1 triangle
            pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
            
            Is = Is(1220:1350,1:140);
            I = I(1220:1350,1:140);
            BW = BW(1220:1350,1:140);
            
            BW = imclearborder([zeros(size(BW,1),1),BW]);
            BW(:,1) = [];
        case 2 % Nuclei clump with 3 nuclei, 1 triangle
            pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
            
            Is = Is(1208:1275,85:170);
            I = I(1208:1275,85:170);
            BW = BW(1208:1275,85:170);
            
            BW = imclearborder(BW);
        case 3 % Nuclei clump with 9 nuclei, 2 connected triangles
            % This is the same as case 1, expect we manually add in another
            % nuclei *to the BW image* to test two connected triangles
            pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
            
            Is = Is(1220:1350,1:140);
            I = I(1220:1350,1:140);
            BW = BW(1220:1350,1:140);
            
            BW = imclearborder([zeros(size(BW,1),1),BW]);
            BW(:,1) = [];
            
            [y,x] = size(BW);
            [Y,X] = meshgrid(1:y,1:x);
            
            
            cent = [30,76];
            r = 20;
            
            R = sqrt((Y-cent(1)).^2 + (X-cent(2)).^2);
            R = R' < r;
            
            BW = BW | R;
        case 4 % Entire nuclei image Cancer
            pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
            
        case 5 % Entire nuclei image OKF
            pth = 'K:\RadiationLab\Cell_Incubation\OKF\OKF_P24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\OKF\Original\OKF_P24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
        case 6 % This is case 1 with the hole filled in
            pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
            if ~exist(pth,'file')
                pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
            end
            I = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
            Is = imfilter(I,fspecial('gaussian',7,1));
            BW = segmentImage(Is);
            
            Is = Is(1220:1350,1:140);
            I = I(1220:1350,1:140);
            BW = BW(1220:1350,1:140);
            
            BW = imclearborder([zeros(size(BW,1),1),BW]);
            BW(:,1) = [];
            BW = imfill(BW,'holes');
        otherwise
            error('Specified case doesn''t exsist')
    end
elseif nargin == 2
    pth = 'K:\Google_Drive\MATLAB\radiationLab\Projects\Cancer_OKF_Cell_Incubation\nucleiClusterImages\';
    
    I = imread([pth,'DAPI', caseNum, '_', cellType, 'P24_image.tif']);
    BW = imread([pth,'DAPI', caseNum, '_', cellType, 'P24_mask.tif']);
    Is = imfilter(I,fspecial('gaussian',7,1));
end

end