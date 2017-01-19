function createResultsData

% Get test image --------------------------------------------------------

[I,~,BW] = getDeclumpTestCase('_2through7','cancer');

% Initialize options ----------------------------------------------------

options = declumpOptions();

options.Max_Radius = 35;
options.Min_Angle = 0.5;

options.Potential_Depth = -1;
options.Potential_Minimum_Location = 2;
options.Potential_Extent = 10;

options.Minimum_Hole_Size = 10;

options.Use_GPU = true;
options.Use_Parallel = true;

options.Debug = true;

selectionTypes = {'curvatureUniformRandom','uniformRandom'};
rs = 5:5:25;
Ns = [150,150];



for selType = 1
    fprintf('%s >> Point selection method = %s...\n', datestr(now,31),selectionTypes{selType})
    options.Point_Selection_Method = selectionTypes{selType};
    N = Ns(selType);
    Info = struct('r0',{},'r_end',{},'centers',{},'solverTime',[],'cutCalculationTime',[],'error',[],'totalComputationTime',NaN);
    Info(N,numel(rs)).totalComputationTime = NaN;
    for ri = 1:numel(rs)
        fprintf('  %s >> WignerSeitzRadius = %g (%d/%d)...\n', datestr(now,31),rs(ri),ri,numel(rs))
        options.Wigner_Seitz_Radius = rs(ri);
        for ni = 1:N
            if ~mod(ni,10)
                fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
            end
            start = tic;
            [~,~,runInfo] = declumpNuclei(I,BW,options);
            totalTime = toc(start);           
            
            Info(ni,ri).r0 = cellfun(@(x) x.r0, runInfo,'UniformOutput',0);
            Info(ni,ri).r_end = cellfun(@(x) x.r_end, runInfo,'UniformOutput',0);
            Info(ni,ri).centers = cellfun(@(x) x.centers, runInfo,'UniformOutput',0);
            Info(ni,ri).solverTime = cellfun(@(x) x.solverTime, runInfo);
            Info(ni,ri).cutCalculationTime = cellfun(@(x) x.cutCalculationTime, runInfo);
            Info(ni,ri).error = cellfun(@(x) x.error, runInfo,'UniformOutput',0);
            Info(ni,ri).totalComputationTime = totalTime;
        end
    end

    pth = 'K:\Google_Drive\MATLAB\radiationLab\Projects\nucleiDeclumpByEM\';
    save([pth, 'results150_rs5by5to25_VintExtent10_' selectionTypes{selType} '_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

end

fprintf('%s >> Done!\n', datestr(now,31))