% createObjectImages ................................................ x
% computeObjectCenters .............................................. x
% computeInitialPoints .............................................. x
% decimateData ................. Fix error when only one edge .......
% interactingParticleSystem ......................................... x
% interactingParticleSystem_convergeEvent ........................... x
% createObjectPartitions ............................................ x
% calculateCuts ..................................................... x
% rebuildBoundary .............. bit messy/not good help ............ x
% orderCutVertices .................................................. x
% computeCutIndices ................................................. x
% computeBoundaryIndices ............................................ x
% computeCutAdjacency ............................................... x
% getPdistInds ...................................................... x
% constrainedObjectCenterTriangulation .............................. x
% findAssociatedCuts ................................................ x
% getSearchInds ..................................................... x
% optimizeCuts ...................................................... x
% checkCutIntegrity ................................................. x
% createTriangleCuts ... need to remove plotting code after test .... x
% chooseWinningCuts ................................................. x



%%
% Get test image --------------------------------------------------------

[I,~,BW] = getDeclumpTestCase('3','Cancer');

% Initialize options ----------------------------------------------------

options = declumpOptions();

options.Max_Radius = 35;
options.Min_Angle = 0.5;

options.Wigner_Seitz_Radius = 20;
options.Potential_Depth = -1;
options.Potential_Minimum_Location = 2;
options.Potential_Extent = 15;

options.Point_Selection_Method = 'curvatureUniformRandom';

options.Use_GPU = true;
options.Use_Parallel = true;
% options.Object_Of_Interest = 234;
options.Debug = true;


start = tic;
[declumpBW,cuts,Info2] = declumpNuclei(I,BW,options);
toc(start)


%%
figure
imshow(declumpBW)

nonEmptyInfo = cellfun(@(x) isfield(x,'centers'),Info2);

centers = cellfun(@(x) x.centers, Info2(nonEmptyInfo),'UniformOutput',false);
centers = cat(1,centers{:});
line(centers(:,2),centers(:,1),'Marker','o','LineStyle','none','Color','b','MarkerSize',4)
% line(cents(:,2),cents(:,3),'Marker','o','MarkerFaceColor','r','LineStyle','none','Color','r','MarkerSize',3)


%%
CC = bwconncomp(BW,8);
L = labelmatrix(CC);

labels = labelCorrectSegmentation(I,declumpedBW,CC);
%%

fig = figure(1);
clf(1)
ax = axes('Parent',fig,'pos',[0 0 1 1]);

% imshow(I,[],'Parent',ax);
cols = {'g','r'};
B = bwboundaries(declumpBW,8,'noholes');
N = numel(B);
bndrySizes = cellfun(@(x) size(x,1), B);


bndry_b = zeros(sum(bndrySizes), 1);
bndry_g = zeros(sum(bndrySizes), 1);
bndry_junk = zeros(sum(bndrySizes), 1);

counter_g = 1;
counter_b = 1;
counter_junk = 1;

for i = 1:N
    try
        linind = B{i}(:,1) + (B{i}(:,2)-1)*size(I,1);
        ind = unique(L(linind));
        ind(ind==0) = [];
%     ind = L(B{i}(1,1),B{i}(1,2));
    if labels(ind,2)==1
        bndry_g(counter_g:bndrySizes(i)+counter_g-1) = B{i}(:,1) + (B{i}(:,2)-1)*size(I,1);
        counter_g = bndrySizes(i)+counter_g;
    elseif labels(ind,2)==3
        bndry_b(counter_b:bndrySizes(i)+counter_b-1) = B{i}(:,1) + (B{i}(:,2)-1)*size(I,1);
        counter_b = bndrySizes(i)+counter_b;
    elseif labels(ind,2)==2
        bndry_junk(counter_junk:bndrySizes(i)+counter_junk-1) = B{i}(:,1) + (B{i}(:,2)-1)*size(I,1);
        counter_junk = bndrySizes(i)+counter_junk;
    end
    catch ME
        i
        [B{i}(1,1),B{i}(1,2)]
        rethrow(ME)
    end
end

bndry_g(bndry_g==0) = [];
bndry_b(bndry_b==0) = [];
bndry_junk(bndry_junk==0) = [];

Irgb = double(I);
Irgb = Irgb-min(Irgb(:));
Irgb = Irgb/max(Irgb(:));
Irgb_g = Irgb;
Irgb_b = Irgb;
Irgb_g(bndry_g) = 1;
Irgb_b(bndry_b) = 1;
Irgb = cat(3,Irgb_b,Irgb_g,Irgb);

imshow(Irgb,'Parent',ax)
%%
ax.Visible = 'off';
pth = 'K:\Google_Drive\NotreDame\RadiationLab\Projects\nucleiDeclumping\Manuscript\Figures\';
imwrite(Irgb,[pth, 'segmentationResults3.png'])
% export_fig(fig,'-native', '-png', [pth, 'segmentationResults2'])

%% Create results data

% Get test image --------------------------------------------------------

[I,~,BW] = getDeclumpTestCase('_2through7','cancer');

% Initialize options ----------------------------------------------------

options = declumpOptions();

options.Max_Radius = 35;
options.Min_Angle = 0.5;

options.Potential_Depth = -1;
options.Potential_Minimum_Location = 2;
options.Potential_Extent = 15;

options.Use_GPU = true;
options.Use_Parallel = true;

options.Debug = true;

rs = 5:5:25;
N = 50;

Info = struct('r0',{},'r_end',{},'centers',{},'cuts',[],'solverTime',[],'cutCalculationTime',[],'totalComputationTime',NaN);
Info(N,numel(rs)).totalComputationTime = NaN;

for ri = 1:numel(rs)
    fprintf('%s >> WignerSeitzRadius = %g (%d/%d)...\n', datestr(now,31),rs(ri),ri,numel(rs))
    options.Wigner_Seitz_Radius = rs(ri);
    for ni = 1:N
        if ~mod(ni,10)
            fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
        end
        start = tic;
        [~,cuts,runInfo] = declumpNuclei(I,BW,options);
        totalTime = toc(start);
        
        Info(ni,ri).r0 = cellfun(@(x) x.r0, runInfo,'UniformOutput',0);
        Info(ni,ri).r_end = cellfun(@(x) x.r_end, runInfo,'UniformOutput',0);
        Info(ni,ri).centers = cellfun(@(x) x.centers, runInfo,'UniformOutput',0);
        Info(ni,ri).cuts = cuts;
        Info(ni,ri).solverTime = cellfun(@(x) x.solverTime, runInfo);
        Info(ni,ri).cutCalculationTime = cellfun(@(x) x.cutCalculationTime, runInfo);
        Info(ni,ri).totalComputationTime = totalTime;
    end
end

fprintf('%s >> Done!\n', datestr(now,31))


