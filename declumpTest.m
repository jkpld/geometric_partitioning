close all

% Get a nuclei mask
% pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_DAPI.tif';
% % pth = 'K:\NotreDame\Projects\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
% % pth = 'K:\Cell_Incubation\Cancer\originalImages\p24_DAPI.tif';
% Id = double(imread(pth,'PixelRegion',{[8193 10240],[8193 10240]}));
% Ids = imfilter(Id,fspecial('gaussian',7,1));
% BW = segmentImage(Ids);
% 
% Idssm = Ids(1220:1350,1:140);
% Idsm = Id(1220:1350,1:140);
% BWsm = BW(1220:1350,1:140);
% 
% % Idssm = Ids(1208:1275,85:170);
% % Idsm = Id(1208:1275,85:170);
% % BWsm = BW(1208:1275,85:170);
% BWsm = imclearborder([zeros(size(BWsm,1),1),BWsm]);
% BWsm(:,1) = [];

[Idsm,Idssm,BWsm] = getDeclumpTestCase(1);

imSize = size(BWsm);

density_Dist = 20;
dxdy = imSize/density_Dist;

x = round(linspace(1,imSize(1),dxdy(1)));
y = round(linspace(1,imSize(2),dxdy(2)));
x2 = round(x(1:end-1) + diff(x)/2);
y2 = round(y(1:end-1) + diff(y)/2);

nx = numel(x);
ny = numel(y);

x = x .* ones(ny,1);
y = y' .* ones(1,nx);
x2 = x2 .* ones(ny-1,1);
y2 = y2' .* ones(1,nx-1);

r0 = [x(:),y(:);x2(:),y2(:)];
r0(r0(:,1)>imSize(2) | r0(:,2)>imSize(1),:) = [];

inds = r0(:,2) + (r0(:,1)-1)*imSize(1);

size(r0)
density = 1/(15*15);
disp here
N = sum(BWsm(:)) * density;
BWsm_e = padarray(BWsm,[1 1],0);
BWsm_e = imerode(BWsm_e,strel('disk',3));
BWsm_e(:,[1,end]) = [];
BWsm_e([1,end],:) = [];

sum(BWsm_e(inds))
good = BWsm_e(inds);
disp here
r0(~good,:) = [];
disp here
r0 = fliplr(r0);

figure
imshow(BWsm,[])
hold on
plot(r0(:,2),r0(:,1),'rx')
% plot(x2(:),y2(:),'bx')

%%
% Get boundary, curvature, and shape markers
kappaSmoothingSigma = 2; % Cancer
markerFilterSize = 5; % Cancer
[B,~,numObjs,bndryTplgy] = bwboundaries(BWsm,8);

[kappa,shapeMarkers,negKappaInds,dB,markersXY] = getCurvatureAndShapeMarkers(B{1},imSize,kappaSmoothingSigma,1,numObjs);
% [pixInds,markers] = createShapeMarkerBasins(B{1},markersXY,imSize,markerFilterSize);

% Get possible starting locations;
[i,j] = find(imerode(BWsm,strel('disk',3)));
validInds = [i,j];

blckSize = 20;

edges_x = 0:blckSize:imSize(2)+blckSize;
edges_y = 0:blckSize:imSize(1)+blckSize;

bins = zeros(size(markersXY));

bins(:,2) = discretize(markersXY(:,2),edges_x);
bins(:,1) = discretize(markersXY(:,1),edges_y);

inds = (1:size(markersXY,1))';

r0_inds = accumarray(bins(all(bins>0,2),:),inds(all(bins>0,2)),[],@(x) x(randi(numel(x),1)));
r0_inds = r0_inds(r0_inds>0);
r0 = markersXY(r0_inds,:);

hold on
plot(r0(:,2),r0(:,1),'bo')
%%

[markersXYunq,~,cnt] = unique(markersXY,'rows');
weight = accumarray(cnt,1,[size(markersXYunq,1),1],[],0);

validInds = markersXYunq;

n = histcounts(cnt,0.5:size(markersXYunq,1)+1);
nunq = unique(n);
cols = myColorMap(numel(nunq),'v1');
% Visualize the nuclei mask ---------------------------------------------
figure(1)
imshow(Idsm,[])
for i = 1:size(markersXYunq,1)
    line(markersXYunq(i,2),markersXYunq(i,1),'Marker','.','Color',cols(n(i)==nunq,:),'MarkerSize',5*n(i),'LineStyle','none')
end
goDark(gcf)
%% Compute the confining potential

% Get potential inside well
V_grid = double(bwdist(~BWsm));
V_grid = imfilter(V_grid,fspecial('gaussian',7,1));
V_grid = 1./V_grid;

% Set potential outside well to strongly increasing with distance. Don't
% set to infinity because, when modeling, if a particle happened to go into
% the Inf region during a time step, then the solution will break.
V_grid_out = (bwdist(BWsm)+1).^5;
V_grid(~BWsm) = V_grid_out(~BWsm);

% Visualize the potential -----------------------------------------------
figure(2)
clf(2)
imagesc(V_grid)
daspect([1 1 1])
clim([0,1])
colorbar

%% Compute the force of the confining potential
[dVx,dVy] = gradient(V_grid);
[dVxx,dVxy] = gradient(dVx);
[~,dVyy] = gradient(dVy);

dV = sqrt(dVx.^2 + dVy.^2);
phase = atan2(-dVy,-dVx);
phase(phase<0) = phase(phase<0)+2*pi;
phase(~BWsm) = 7;
% Visualize the force ---------------------------------------------------
figure(3)
clf(3)
imagesc(dV)
clim([0,1])
daspect([1 1 1])
colorbar

% Visualize the direction of the force ----------------------------------
% First create a cyclical colormap

cmap = brewermap(36,'GnBu');
cmapP = [cmap; flipud(cmap);0 0 0];

% Now plot the direction
figure(4)
clf(4)
imagesc(phase)
daspect([1 1 1])
colormap(cmapP)
colorbar

% Create interpolating functions for confining force and potential

[Y,X] = ndgrid(1:size(BWsm,1),1:size(BWsm,2));

dVx = griddedInterpolant(Y,X,dVx);
dVy = griddedInterpolant(Y,X,dVy);

dVxx = griddedInterpolant(Y,X,dVxx);
dVxy = griddedInterpolant(Y,X,dVxy);
dVyy = griddedInterpolant(Y,X,dVyy);

V = griddedInterpolant(Y,X,V_grid);

% dVx = @(x) interp2(X,Y,dVx,x(:,2),x(:,1),'*linear');
% dVy = @(x) interp2(X,Y,dVy,x(:,2),x(:,1),'*linear');

%% Create particle interaction force function

% Short range attraction long range repulsion ---------------------
% V_LJ = @(r,ep,sg,a) 4*ep*( (sg./r).^(2*a) - (sg./r).^a );
% dV_LJ = @(r,ep,sg,a) -4*a*ep*( 2*(sg./r).^(2*a) - (sg./r).^a ).*(1./r);
%
% ep = 10;
% sg = 1;
% a = 2;
%
% Vint = @(D) V_LJ(D,ep,sg,2*a) - V_LJ(D,0.02*ep,2*sg,a/2);
% dVint = @(D) dV_LJ(D,ep,sg,2*a) - dV_LJ(D,0.02*ep,2*sg,a/2);

A = 1;
sig = 3;
mu = 2;

Vint = @(D) 1./D - A * exp(-(D-mu).^2/(2*sig^2));
dVint = @(D) -1./D.^2 + (A*(D-mu)/(sig^2)) .* exp(-(D-mu).^2/(2*sig^2));
% dVint = @(D) -1./D.^2 + (1*(D-2)/(16)) .* exp(-(D-2).^2/(32));

% Coulomb ---------------------------------------------------------
dVint = @(D) -1./D.^2; % Coulomb force
Vint = @(D) 1./D; % Coulomb potential


%% Initialize particles

% Set the number of particles ==========================================
% N = 50;%round(size(validInds,1)/7)
N = size(r0,1);
% N = 8;
% Package the needed inputs ============================================

q = N^(-1/3);%ones(N,1)*N^(-1/3);
m = 1;%ones(N,1);%@(t) ones(N,1)+0.05*(t(:))';%ones(N,1);
alpha = @(t) zeros(N,1) + t*5e-4; %2e-4 for coulomb

extraInputs.q = q;
extraInputs.m = m;
extraInputs.alpha = alpha;
extraInputs.dVx = dVy; % Note that the order of dVx dVy is switched. This is because y is the first dimension and x is the second.
extraInputs.dVy = dVx;
extraInputs.dVint = dVint;

% Get the pdist indices needed
N = 5;
pdistInds = getPdistInds(N);
accumIndsOut = sparse(1:N).' == pdistInds(:,1).'; % sparse matrix with a full size of N x N*(N-1)/2
accumIndsIn  = sparse(1:N).' == pdistInds(:,2).'; % sparse matrix with a full size of N x N*(N-1)/2

N_pair = N*(N-1)/2;

S = sparse(1:N_pair,pdistInds(:,1),1,N_pair,N,N_pair);
S2 = sparse(1:N_pair,pdistInds(:,2),1,N_pair,N,N_pair);

createRpairs = S - S2;
% accumIn = S2';
% accumOut = S';
% accumIn - accumeOut = (S2-S)' = -creatRpairs'

extraInputs.pdistInds = pdistInds;
extraInputs.accumIndsOut = accumIndsOut;
extraInputs.accumIndsIn = accumIndsIn;
extraInputs.accumInds = accumIndsIn-accumIndsOut;
extraInputs.createRpairs = createRpairs;

extraInputsE.q = q;
extraInputsE.m = m;
extraInputsE.alpha = alpha;
extraInputsE.V = V;
extraInputsE.Vint = Vint;
extraInputsE.pdistInds = pdistInds;
extraInputsE.accumIndsOut = accumIndsOut;
extraInputsE.accumIndsIn = accumIndsIn;
extraInputsE.accumInds = accumIndsIn-accumIndsOut;

extraInputsJ.q = q;
extraInputsJ.m = m;
extraInputsJ.alpha = alpha;
extraInputsJ.dVxx = dVyy; % Note that the order of dVx dVy is switched. This is because y is the first dimension and x is the second.
extraInputsJ.dVxy = dVxy;
extraInputsJ.dVyy = dVxx;
extraInputsJ.pdistInds = pdistInds;

extraInputsM.m = @(t) ones(N,1)+0.01*t;
extraInputsM.N = N;
% Create initial conditions ============================================
% w0 = weight;
% w0 = 1./weight;
% w0 = ones(size(weight));
w0 = ones(size(validInds,1));
w0 = w0/sum(w0);
% r0 = datasample(validInds,N,1,'Weights',w0,'Replace',false);
v0 = 0.01*ones(N,2);

offset = (0:4:N*4-1)';

rInds = [1,2] + offset;
vInds = [3,4] + offset;

y0 = zeros(N,1);
y0(rInds) = r0;
y0(vInds) = v0;


% J = interactingParticleSystemJacobian(0,y0,extraInputsJ);

% Get the energy of the initial conditions.
E = electronSystemEnergy([],y0,V,Vint,pdistInds);
E0 = sum(E);
fprintf('Initial Energy = %0.3g\n',E0)

% Set up the ODE solver ================================================

odeFun = @(t,y) interactingParticleSystem(t,y,extraInputs);
jacobianFun = @(t,y) interactingParticleSystemJacobian(t,y,extraInputsJ);
massFun = @(t) interactingParticleSystemMass(t,extraInputsM);
convergeFun = @(t,y) electronSystemConverge_event(t,y);
options = odeset('Events',convergeFun,...
    'Vectorized','on',...
    'RelTol',1e-4,...
    'AbsTol',1e-6,...
    'NormControl','on',...
    'Jacobian',jacobianFun);%,...
%     'Mass',massFun,...
%     'MStateDependence','none');

% odeFun = @(t,y) electronSystem(t,y,dVy,dVx,dVint,pdistInds); % Note that the order of dVx dVy is switched. This is because y is the first dimension and x is the second.
% options = odeset('Events',@electronSystemConverge_event);

% Solve ODE ============================================================

% NumKicks = 1;
%
% Times = [];
% Sols = [];
%
% for i = 1:NumKicks
%     [T,Y] = ode23s(odeFun,0:10:200,y0,options);
%
%     % Give the system a little kick.
%     y0_2 = Y(end,:)';
%     y0_2_v = y0_2(vInds);
%     y0_2_v = y0_2_v + 0.5*(rand(N,2)-0.5);
%     y0_2(vInds) = y0_2_v;
% end

start = tic;


% Solve the system for a short time
profile on
tic
[T,Y,te,ye,ie] = ode23(odeFun,0:10:1500,y0,options);
% [T,Y,te,ye,ie] = ode23(odeFun,[0,1500],y0,options);
toc

% % Give the system a little kick.
% P = sum(interactingParticleSystemPotentialEnergy([],Y,extraInputsE),2);
% [~,resetInd] = min(P);
% y0_2 = Y(resetInd,:)';
% % y0_2_v = y0_2(vInds);
% ur = rand(N,2)-0.5;
% ur = ur ./ sqrt(sum(ur.^2,2));
% y0_2_v = 0.2*ur;%2*sqrt(sum(y0_2_v.^2,2)) .* ur;%0.8*y0_2_v + 0.2*ur;%
% y0_2(vInds) = y0_2_v;
%
% tic
% [T2,Y2,te,ye,ie] = ode45(odeFun,0:10:300,y0_2,options);
% toc
%
% % Give the system another little kick.
% P = sum(interactingParticleSystemPotentialEnergy([],Y2,extraInputsE),2);
% [~,resetInd] = min(P);
% y0_2 = Y2(resetInd,:)';
% % y0_2_v = y0_2(vInds);
% ur = rand(N,2)-0.5;
% ur = ur ./ sqrt(sum(ur.^2,2));
% y0_2_v = 0.1*sqrt(sum(y0_2_v.^2,2)) .* ur;%0.9*y0_2_v + 0.1*ur;%0.3*ur;%
% y0_2(vInds) = y0_2_v;
%
% % Solve the system out
% tic
% [T3,Y3,te,ye,ie] = ode45(odeFun,0:10:1000,y0_2,options);
% toc
toc(start)

profile off
profile viewer


fprintf('Final Energy = %0.3g\n', sum(electronSystemEnergy([],Y(end,:)',V,Vint,pdistInds)))

% T = [T;T2+200;T3+500];
% Y = [Y;Y2;Y3];

% Get the total energy as function of time
E = interactingParticleSystemEnergy(T,Y,extraInputsE);
P = interactingParticleSystemPotentialEnergy([],Y,extraInputsE);
P = sum(P,2);
E = sum(E,2);

% Visualize particle movement ===========================================

try close(fig5), catch, end;
fig5 = figure(5);
ax5 = axes('Parent',fig5);
contourf(log(V_grid),linspace(-3,3,20),'Parent',ax5)
set(gca,'ydir','reverse')
cmapBW = createColorGradient([1 1 1],[0 0 0],56);
colormap(cmapBW)
% clim([0,0.5])
colorbar

particleColors = myColorMap(N,'v1');

clear l
for i = N:-1:1
    l(i) = line(NaN,NaN,1,'Marker','o','Color',particleColors(i,:),'MarkerFaceColor',particleColors(i,:),'LineStyle','none','Parent',ax5);
end

for i = 1:numel(T)
    for j = 1:N
        l(j).YData = Y(i,rInds(j,1));
        l(j).XData = Y(i,rInds(j,2));
    end
    pause(0.05)
end
fprintf('finished!\n')

pause(0.5)

% Visualize particle speeds
try close(fig6), catch, end
fig6 = figure(6);
ax6 = axes('Parent',fig6);

for j = 1:N
    line(T,sqrt(sum(Y(:,vInds(j,:)).^2,2)),'Color',particleColors(j,:),'LineWidth',2,'Parent',ax6)
end
xlabel('Time')
ylabel('Particle speed')
drawnow;
goDark(gcf)

% Visualize total and potential energy
try close(fig7), catch, end
fig7 = figure(7);
ax7 = axes('Parent',fig7);

line(T,E,'Color','b','LineWidth',2,'Parent',ax7)
line(T,P,'Color','r','LineWidth',2,'Parent',ax7)
xlabel('Time')
ylabel('Energy')
legend({'Total energy','Potential energy'})
drawnow;
goDark(gcf)
figure(5)

% Visualize mean particle speed
% try close(fig8), catch, end
% fig8 = figure(8);
% ax8 = axes('Parent',fig8);
% speed = zeros(size(T,1),1);
% for j = 1:N
%     speed = speed + sqrt(sum(Y(:,vInds(j,:)).^2,2));
% end
% line(T,speed/N,'Color',particleColors(j,:),'LineWidth',2,'Parent',ax8)
% xlabel('Time')
% ylabel('Mean particle speed')
% drawnow;
% goDark(gcf)

%% Form mask

y_end = Y(end,:)';
r_cents = round(y_end(rInds));

r_centsInd = r_cents(:,1) + (r_cents(:,2)-1) * size(BWsm,1);

mask = zeros(size(BWsm));
mask(r_centsInd) = 1;
mask = imdilate(mask,strel('disk',2));
% mask = imerode(mask,strel('disk',1));
CC = bwconncomp(mask);
props = regionprops(CC,'Centroid');
cents = fliplr(cat(1,props.Centroid));
% %%
try close(fig9), catch, end
fig9 = figure(9);
imshow(mask)


D = -bwdist(~BWsm);
D = imimposemin(D,mask);
% D = bwdist(mask);
% D(~BWsm) = inf;

try close(fig10), catch, end
fig10 = figure(10);
imagesc(D)
daspect([1 1 1])
colorbar
drawnow;
goDark(gcf)


L = watershed(D,8);
BWmask = (L>0) & BWsm;


try close(fig11), catch, end
fig11 = figure(11);
imagesc(BWmask)
% daspect([1 1 1])
% colorbar
% drawnow;
% goDark(gcf)

%%

bndryB = cell(1,numObjs);
dbndry = cell(1,numObjs);
markersXY = cell(1,numObjs);
kappa = cell(1,numObjs);


for i = 1:numObjs
    
    Prnt = B{i};
    if ~ispolycw(B{i}(:,1),B{i}(:,2))
        Prnt = flip(Prnt,1);
    end
    
    [kappaP,~,~,dBP,markersXYP] = getCurvatureAndShapeMarkers(Prnt,imSize,kappaSmoothingSigma,[],[]);
    Prnt(end,:) = [];%nan;
    
    children = bndryTplgy(:,i);
    if any(children)
        children = find(children);
        numChildren = numel(children);
        
        child = cell(1,numChildren);
        kappa_child = cell(1,numChildren);
        dB_child = cell(1,numChildren);
        markersXY_child = cell(1,numChildren);
        
        for j = 1:numChildren
            
            child{j} = B{children(j)};
            if ispolycw(child{j}(:,1),child{j}(:,2))
                child{j} = flip(child{j});
            end
            
            [kappa_child{j},~,~,dB_child{j},markersXY_child{j}] = getCurvatureAndShapeMarkers(child{j},imSize,kappaSmoothingSigma,[],[]);
            child{j}(end,:) = nan;
        end
        
        bndryB{i} = [Prnt; nan, nan; cat(1,child{:})];
        kappa{i} = [kappaP; cat(1,kappa_child{:})];
        dbndry{i} = [dBP; cat(1,dB_child{:})];
        markersXY = [markersXYP; cat(1,markersXY_child{:})];
        
    else
        bndryB{i} = [Prnt; nan, nan];
        kappa{i} = kappaP;
        dbndry{i} = dBP;
        markersXY = markersXYP;
    end
    
    % Remove the row of NaN's at the end.
    bndryB{i}(end,:) = [];
    kappa{i}(end,:) = [];
    dbndry{i}(end,:) = [];
    
end



%% Compute declumping edges

options.maxRadius = 30;
options.minAngle = 0.5;
options.searchRadius = 15;

bndry = bndryB{1};
dB = dbndry{1};

% holdInd = find(isnan(bndry(:,1)));
% bndry(holdInd:end,:) = [];
% dB(holdInd:end,:) = [];

% Get boundary normals (pointing inward)
n = -[-dB(:,2), dB(:,1)]; % Mx2
n = n ./ sqrt(sum(n.^2,2)); % Mx2


cnts = cents; % Nx2
% cnts(5,:) = [];
tic
% Create initial object partitions and cuts
[cuts,ncuts,cutCntrs,partitions,Info] = createObjectPartitions(bndry,n,cnts,options);

assignedCenter = Info.assignedCenter;

% Split up the boundary into a cell array of seperate contours and get the
% cut vertex indices in each seperate contour
[B_cell,n_cell,cutVertices_cell,cutVertices,cutVerticesContourNumber,B_contourNumber] = ...
    splitContoursAndGetIndices(bndry,n,cuts,ncuts);

% Get the unpaired cuts.
unPaired = isnan(cutCntrs(:,2));
unPairedPartitionNumber = cutCntrs(unPaired,1);

% Create the unpaired edges data for optimization with optimizeUnpairedCuts
unpairedCuts = [cutVertices(unPaired,:),cutVerticesContourNumber(unPaired,:)];

% Optimize the unpaired cuts
[optimizedEdges,optimizedEdgeNorms,cutVertices_cell] = optimizeUnpairedCuts(unpairedCuts,B_cell,n_cell,cutVertices_cell,options);

cuts(unPaired,:) = optimizedEdges;
ncuts(unPaired,:) = optimizedEdgeNorms;
toc

figure(101)
clf(101)
line(bndry(:,2),bndry(:,1))
hold on
quiver(bndry(:,2),bndry(:,1),n(:,2),n(:,1))
set(gca,'ydir','reverse')

for i = 1:size(bndry,1)
    if ~isnan(assignedCenter(i))
        line([bndry(i,2),cnts(assignedCenter(i),2)], [bndry(i,1),cnts(assignedCenter(i),1)],'color',0.6*[1 1 1])
    end
end

for i = 1:size(cnts,1)
    line(partitions{i}(:,2),partitions{i}(:,1),'Color','g')
end

for i = 1:size(cuts,1)
    line(cuts(i,[2,4]),cuts(i,[1,3]),'color','r')
end

plot(cnts(:,2),cnts(:,1),'rx')
daspect([1 1 1])
drawnow;
goDark(gcf)

for i = 1:size(cnts,1)
    text(cnts(i,2),cnts(i,1),num2str(i),'FontSize',15,'Color',[0 1 0],'FontWeight','bold')
end

%%
%{
  This is a comment! That's awesome!
%}

%

triangleList = delaunay(cnts);

triEdges = [triangleList(:,[1,2]); triangleList(:,[2,3]); triangleList(:,[3,1])];
triEdges = unique(sort(triEdges, 2, 'ascend'),'rows');

xVerts = [cnts(triEdges(:,1),1) cnts(triEdges(:,2),1) nan(size(triEdges,1),1)]';
yVerts = [cnts(triEdges(:,1),2) cnts(triEdges(:,2),2) nan(size(triEdges,1),1)]';

[~,~,intrsctns] = intersections(xVerts(:),yVerts(:),bndry(:,1),bndry(:,2));

badEdges = triEdges(unique(ceil(intrsctns/3)),:);

for i = 1:size(badEdges,1)
    triangleList(any(triangleList == badEdges(i,1),2) & any(triangleList == badEdges(i,2),2),:) = [];
end

% For each triangle left over
% * Get the center of the triangle
% * Get the center of the three edges
% * Generate three new edges that go from the triangle center, through the
% edge centers, all the way to the boundary.
% * Determine what two edges we are trying to improve
%  # The row of the triangle list will contain three center indices.
%  # Find the two edges in "cuts" whos have combinations of the indices in
%  the row of the triangle list.
% * If any of the new edges intersects with an edge that we are not trying
% to replace, then remove this triangle.
% * If the new edges are all still valid, then compute mean overlap of the
% new edges with the image edge (weighed average of the first and second
% derivatives of the image intensity), and compute the mean overlap of the
% old edges with the image edge.
% * Use the set of edges that have the largest mean overlap.

% Get the maximum distance a point could be from the boundary.
MAX_CUT_LENGTH = sqrt(sum(range(bndry).^2));

% Get the cut line segments
cutLines = reshape([cuts, nan(size(cuts,1),2)]',2,size(cuts,1)*3)';

% Create the set of pixels that the cuts cover.
cutInds = cell(size(cuts,1),1);
for i = 1:size(cuts,1)
    inds = rayTrace(cuts(i,1:2),cuts(i,3:4));
    cutInds{i} = inds(:,1) + (inds(:,2)-1)*imSize(1);
end

% Create a search range index array;
searchRange = -10:10;

% Get image edges
S = getImageEdges(Idssm);

% Find centers that have more than one unpaired edge
unpairedCuts = cutCntrs(isnan(cutCntrs(:,2)),1);
uniqueUnpairedEdges = unique(unpairedCuts);
binEdges = [uniqueUnpairedEdges(1)-0.5; uniqueUnpairedEdges(1:end-1) + diff(uniqueUnpairedEdges)/2; uniqueUnpairedEdges(end)+0.5];
unpairedEdgesCount = histcountsmex(unpairedCuts,binEdges);
centersWithMoreThanOneUnpairedEdge = uniqueUnpairedEdges(unpairedEdgesCount>1);

% Get the length of each boundary contour
B_cell_lengths = cellfun(@(x) size(x,1),B_cell);
offsets_boundary_to_cell = cumsum(B_cell_lengths+1);
offsets_boundary_to_cell = [0; offsets_boundary_to_cell(1:end)];

% Create an array of binEdges for use later
binEdges = 0.5:3.5;

newCuts = cell(size(triangleList,1)*3,1);
oldCutsToRemove = false(size(cnts,1),1);

debug = true;
counter = 1;
if ~isempty(triangleList)
    for i = size(triangleList,1):-1:1
        i
        % ===============================================================
        % Part 1:
        % * Get center of triangle
        % * Get mid point of each edge of the triangle               
        % * Generate three new proposed cuts that go from the triangle
        %   center, through the edge mid points, all the way to the
        %   boundary.
        % * Intersects the three new proposed cuts with all cuts previously
        %   found.
        
        % Get the center of the triangle
        triCent = mean(cnts(triangleList(i,:),:),1);
        
        % Get the centers of the triangle edges
        edgCents = mean(cat(3,cnts(triangleList(i,:),:),cnts(circshift(triangleList(i,:),1),:)),3);
        
        % Get the unit vectors that start at the center and go through the
        % triangle's edge centers.
        uv = edgCents-triCent;
        uv = uv ./ sqrt(sum(uv.^2,2));
        
        % Get the end points of the new propopsed edges.
        edgEnds = triCent + uv*MAX_CUT_LENGTH;
        
        % Align the three new proposed edges as a NaN delimeted array.
        xTmp = [triCent(1)*[1;1;1], edgEnds(:,1) nan(3,1)]';
        yTmp = [triCent(2)*[1;1;1], edgEnds(:,2) nan(3,1)]';
        
        % Get the intersections between the new proposed edges and the old
        % set of edges.
        [~,~,intrsctns] = intersections(cutLines(:,1),cutLines(:,2),xTmp(:),yTmp(:));
        intersectedEdgeInds = unique(ceil(intrsctns/3));
        % End Part 1 ====================================================
        
        % ===============================================================
        % Part 2:
        %         
        % Need to find all cuts associated with this triangle. An
        % associated cut could be a "hard edge" (only one cut seperateing
        % the two centers), a "soft edge" (two cuts seperating two
        % centers), or an "unpaired edge" that is only associated with one
        % center. If there is a soft or hard edge between two centers then
        % take that edge, if there is not a soft or hard edge associated
        % with one of the centers, then there must be an unpaired edge, and
        % we need to find the unpaired edge closest to the triangle center.
        
        % Get paired (hard and soft) edges.
        currentCuts = sum(sum(cutCntrs == permute(triangleList(i,:),[1,3,2]),3),2)==2;
        
        % If this set of currentEdges does not contain reference to one of
        % the triangle centers, then look for unpaired edges.
        cntrsNotIncluded = ~ismember(triangleList(i,:),unique(cutCntrs(currentCuts,:)));
        
        if any(cntrsNotIncluded)
            % Note that unpaired edges are always listed in the first
            % column of cutCtnrs.
            unpairedCuts = any(cutCntrs(:,1) == triangleList(i,cntrsNotIncluded),2) & isnan(cutCntrs(:,2));

            % Check to see if a center of the triangle has more than one
            % unpaired edges, and if it does, only use the edge whose midpoint
            % is closest to the triangle center.
            moreThanOneUnpaired = find(any(triangleList(i,:) == centersWithMoreThanOneUnpairedEdge,2)); %nx3 -> nx1
            
            for j = 1:numel(moreThanOneUnpaired)
                
                tmpCntr = centersWithMoreThanOneUnpairedEdge(moreThanOneUnpaired(j));
                tmpEdgeInds = find(cutCntrs(:,1) == tmpCntr);
                tmpEdges = cuts(tmpEdgeInds,:);
                tmpEdgesMidPoints = (tmpEdges(:,1:2) + tmpEdges(:,3:4))/2;
                
                [~,edgeToUse] = min( sum((tmpEdgesMidPoints - triCent).^2,2) );
                tmpEdgeInds(edgeToUse) = [];
                unpairedCuts(tmpEdgeInds) = false;
            end
        else
            unpairedCuts = false(size(currentCuts));
        end
        
        % Combine the cuts to get the complete set.
        currentCuts = currentCuts | unpairedCuts;
        % End Part 2 ====================================================
        

        
        % ===============================================================
        % Part 3:
        %
        % * Ensure that if the new propsed edges intersect with the
        %   previous edges, that they only intersect with the two edges we
        %   are proposing to replace. If the new proposed edges intersect
        %   with other edges, then do not consider these new proposed
        %   edges.
        % * Optimize the new proposed cut's intersection with the boundary.
        %   # First remove the currentVertices from consideration.
        %   # Use all other (not current) vertices to bound the possible
        %     motion of the new porposed cut's intersections.
        
        if ~any(~ismember(intersectedEdgeInds, find(currentCuts)))
            
            % Find closest intersection of new proposed cuts with boundary.
            [~,~,intrsctns,intrsctns2] = intersections(xTmp(:),yTmp(:),bndry(:,1),bndry(:,2)); % !! first offset the boundary by the normal vector by 0.5 !!
            
            if numel(intrsctns) > 3
                % If there are more than three intersectinos then we need
                % to use the intersections from the proposed cuts.
                
                % Find the intersection of each proposed cut that is
                % closest to the triangle center.
                proposedCutIntersections = ceil(intrsctns/3);
                closestIntersection = accumarray(proposedCutIntersections, intrsctns, [3,1], @min, nan);
                if any(isnan(closestIntersection))
                    error('A NaN was returned when finding the closest intersection between the proposed triangle cuts and the boundary.')
                end
                % Compute the fractional distance along each proposed cut
                % to the closest intersection.
                closestIntersection = closestIntersection - floor(closestIntersection);

                % Find the xy location of the closest intersections.
                closestIntersectionLocations = triCent + (MAX_CUT_LENGTH * closestIntersection) .* uv;
                
                % Find the closest boundary vertex to these intersection
                % locations.
                [~,newBinds] = min(sum((bndry - permute(closestIntersectionLocations,[3,2,1])).^2,2));
                newBinds = permute(newBinds,[3,2,1]);
            else
                % If there are only three intersections, then we can get
                % the intersection location directly from the boundary
                % intersections.
                newBinds = round(intrsctns2);
            end
            
            if debug
                TriangleInfo(i).originalCutVertices = bndry(newBinds,:);
            end
            
            % Now that we have the bounary indices of the proposed cuts
            % intersectios with the boundary, we need to optimize the
            % positions. 
            
            % * Convert these intersection indices into the indices for the
            %   boundary contours.
            % * Remove the cut vertices of the edges we are trying to
            %   replace from cutVertices and add in the proposed cut edges
            % * Get the indices we will search for an optimized value over
            % * Optimize the intersection by maximizing the dot product of
            %   the boundary normal and the proposed cut direction and by
            %   minimizing the proposed cut distance.

            intersectionContours = B_contourNumber(newBinds);
            newBinds = newBinds - offsets_boundary_to_cell(intersectionContours);
            
            tmpCutVertices = cutVertices(~currentCuts,:); % THIS IS WRONG! CUTVERTICES does not have the new vertices after optimization of unpaired edges.
            tmpCutVerticesContourNumber = cutVerticesContourNumber(~currentCuts,:);
            
            v = accumarray([tmpCutVerticesContourNumber(:); intersectionContours(:)], [tmpCutVertices(:); newBinds(:)], [numel(B_cell),1], @(x) {sort(x)},[]);
            for j = 1:3
                
                contourNumber = intersectionContours(j);
                searchInds = getSearchInds(newBinds(j),B_cell_lengths(contourNumber),v{contourNumber},options.searchRadius);
                
                r1 = B_cell{contourNumber}(searchInds,:);
                n1 = n_cell{contourNumber}(searchInds,:);
                
                r1 = triCent - r1;
                d1 = sqrt(sum(r1.^2,2));
                r1 = r1 ./ d1;
                
                a1 = sum( r1 .* n1, 2);
                
                valueToOptimize = a1 ./ d1;
                [~,bestInd] = max(valueToOptimize);
                newBinds(j) = searchInds(bestInd);
            end
            newBinds = newBinds + offsets_boundary_to_cell(intersectionContours);
            newCutVerts = bndry(newBinds,:);
            newCutVertNorms = n(newBinds,:);

            if debug
                TriangleInfo(i).optimizedCutVertices = newCutVerts;
            end
            
            % End part 3 ================================================
            
            
            % ===========================================================
            % Part 4: Optimize triangle center

            triCent = round(triCent);
            
            if debug
                TriangleInfo(i).originalTriangleCenter = triCent;
            end
            tic
            % ---------------------
            % Run optimization here
            
            % Optimize the center of the triangle over the points inside of
            % a second triangle whose vertices are a fraction along the
            % exsisting cuts.
            TRIANGLE_CENTER_SEARCH_FRAC = 0.7;
            triangleSearchBounds = (1-TRIANGLE_CENTER_SEARCH_FRAC)*triCent + TRIANGLE_CENTER_SEARCH_FRAC * newCutVerts;
            minSearchBounds = round(min(triangleSearchBounds) - 1);
            maxSearchBounds = round(max(triangleSearchBounds) + 1);
            
            % Create the x and y values to search over such that the
            % triangle center is always one of the search values.
            x = [triCent(1):-2:minSearchBounds(:,1), triCent(1)+2:2:maxSearchBounds(:,1)];
            y = [triCent(2):-2:minSearchBounds(:,2), triCent(2)+2:2:maxSearchBounds(:,2)];
            X = x' .* ones(1,numel(y));
            Y = y  .* ones(numel(x),1);
            
            search_r = [X(:),Y(:)];
            
            % Remove any search point that is outside of the search
            % triangle.
            in = inpolygon(search_r(:,1),search_r(:,2),triangleSearchBounds(:,1),triangleSearchBounds(:,2));
            search_r(~in,:) = [];
            
            
            r1 = search_r - permute(newCutVerts,[3,2,1]);
            d1 = sqrt(sum(r1.^2,2));
            r1 = r1 ./ d1;
            
            a1 = sum(r1 .* permute(newCutVertNorms,[3,2,1]),2);
            % Optimize the mean dot product divide by the mean distance. 
            % Note : if you try to optimize the mean of the dot prodcut
            % divided by the distance instead you will not get good results
            % because one let of the triangle will get very short.
            valueToOptimize = sum(a1,3)./sum(d1,3);

            [~,bestTriCentInd] = max(valueToOptimize);
            triCent = search_r(bestTriCentInd,:);

            % ---------------------
            
            if debug
                TriangleInfo(i).triangleCenterSearchPoints = search_r;
                TriangleInfo(i).optimizedTriangleCenter = triCent;
            end
            
            
            % End Part 4 ================================================
            
            
            % ===========================================================
            % Part 5: 
            % Choose between new proposed cuts or original cuts by looking
            % at the cuts overlap with the intensity image's edges.
            
            % Get the linear indices of the new edges in the image
            proposedEdgeLinInds = cell(3,1);
            for j = 1:3
                inds = rayTrace(triCent,newCutVerts(j,:));
                proposedEdgeLinInds{j} = inds(:,1) + (inds(:,2)-1)*imSize(1);
            end
            
            % Find overlap between new and old edges with the image edge.
            overlapOld = mean(S(cat(1,cutInds{currentCuts})));
            overlapNew = mean(S(cat(1,proposedEdgeLinInds{:})));
                        [overlapOld, overlapNew]
                        
             % Select the winning set of edges.
            if overlapNew > overlapOld
                newCuts((1:3) + (counter-1)*3) = proposedEdgeLinInds;
                oldCutsToRemove(currentCuts) = true;
                counter = counter + 1;
            end
            % End Part 5 ================================================
            
        end % If new proposed cuts do not overlap with any cuts they are not supposed to
    end % for each triangle
end % if there are triangles

% Remove any cuts that are to-be-removed and add in any new cuts;
cutInds(oldCutsToRemove) = [];
cutInds = [cutInds; newCuts];

if debug
    for i = 1:numel(TriangleInfo)
        if ~isempty(TriangleInfo(i).optimizedCutVertices)
            for j = 1:3
                
                x = [TriangleInfo(i).originalTriangleCenter(2); TriangleInfo(i).optimizedCutVertices(j,2)];
                y = [TriangleInfo(i).originalTriangleCenter(1); TriangleInfo(i).optimizedCutVertices(j,1)];
                line(x,y,'Color','m')

                x = [TriangleInfo(i).originalTriangleCenter(2); TriangleInfo(i).originalCutVertices(j,2)];
                y = [TriangleInfo(i).originalTriangleCenter(1); TriangleInfo(i).originalCutVertices(j,1)];
                line(x,y,'Color','y')
                
                x = [TriangleInfo(i).optimizedTriangleCenter(2); TriangleInfo(i).optimizedCutVertices(j,2)];
                y = [TriangleInfo(i).optimizedTriangleCenter(1); TriangleInfo(i).optimizedCutVertices(j,1)];
                line(x,y,'Color','c')
            end
            
            plot(TriangleInfo(i).triangleCenterSearchPoints(:,2),TriangleInfo(i).triangleCenterSearchPoints(:,1),'ys','MarkerSize',4);
        end
    end
end

% plot(triCent(:,2),triCent(:,1),'yo')
% plot(newEdges(:,2),newEdges(:,1),'yo')


% BWpart = BWsm;
BWpart = imfill(BWsm,'holes');
BWpart(cat(1,cutInds{:})) = 0;
BWpart = bwareaopen(BWpart,45);
BWpart_b = bwboundaries(BWpart);
figure(102)
imshow(Idsm,[])
for i = 1:numel(BWpart_b)
    line(BWpart_b{i}(:,2),BWpart_b{i}(:,1),'color','g')
end
drawnow;
goDark(gcf)


% %% Just minimize with fminunc
%
% % N = 6;
%
% r0 = validInds(randperm(size(validInds,1),N),:)';
% r0 = r0(:);
% minFun = @(y) electronSystemEnergy_onlyPos(y,V,dVy,dVx,Vint,dVint,pdistInds);
%
% % [E,dE] = minFun(r0)
% %
% options = optimoptions('fminunc','SpecifyObjectiveGradient',false,'CheckGrad',false);
% tic
% r = fminunc(minFun,r0,options);
% % r = simulannealbnd(minFun,r0);
% toc
% r = reshape(r,2,N)';
%
% try close(fig7), catch, end;
% fig7 = figure(7);
% ax7 = axes('Parent',fig7);
% contourf(V_grid,logspace(-2,0,20),'Parent',ax7)
% set(gca,'ydir','reverse')
% cmapBW = createColorGradient([1 1 1],[0 0 0],56);
% colormap(cmapBW)
% colorbar
%
% for i = N:-1:1
%     line(r(i,2),r(i,1),1,'Marker','o','Color',particleColors(i,:),'MarkerFaceColor',particleColors(i,:),'LineStyle','none','Parent',ax7);
% end
%


%%

S2 = abs(imfilter(Idssm,fspecial('log',14,2),'symmetric'));
[Gx,Gy] = smoothGradient(Idssm,1);
S1 = sqrt(Gx.^2 + Gy.^2);
S = S2/prctile(S2(:),80) + S1/prctile(S1(:),80);
S = imclose(S,strel('disk',3));
% S = imopen(S,strel('disk',3));
% S = imclose(S,strel('disk',3));
figure
imshow(S,[])


%% play with potentials

r = 1e-2:0.01:20;

epsilon = 1;
sigma = 1;
alpha = 2;

V_LJ = @(r,epsilon,sigma,alpha) 4*epsilon*( (sigma./r).^(2*alpha) - (sigma./r).^alpha );
Yuk = @(r,A,xi) A*exp(-r/xi)./(r/xi);
V_exp = @(r,mu,sig,A) 1./r - A*exp(-(r-mu).^2/(2*sig^2));
% try close(fig100), catch, end
fig100 = figure(100);
clf(fig100);
% SRLR1 = Yuk(r,epsilon*0.05,sigma*2)+V_LJ(r,epsilon,sigma,alpha);
% SRLR2 = V_LJ(r,epsilon,sigma,alpha*2) - V_LJ(r,0.02*epsilon,sigma*2,alpha/2);
SRLR2 = V_exp(r,2,2,0.7);
SRLR3 = V_exp(r,2,2,1/2);
% line(r,SRLR1,'Color','r','LineWidth',2)
line(r,SRLR2,'Color','g','LineWidth',2)
line(r,SRLR3,'Color','y','LineWidth',2)
line(r,1./r,'Color','b','LineWidth',2)
legend({'1/r - gauss(2,2,0.7)','1/r - gauss(2,2,1+1/mu)','1/r'})

ylim([min([min(SRLR2),min(SRLR3)]),3])
drawnow;
line([0,20],[0 0],'Color',0.5*[1 1 1],'LineStyle','--')
goDark(gcf)

%%

depths = -0.1:-0.1:-2;
centers = 1.5:0.5:4;
extents = 8:1:20;

potentialParameters = zeros(numel(depths),numel(centers),numel(extents),3);

options = optimset('display','none','TolX',1e-7,'TolFun',1e-7,'MaxFunEvals',500);
tic
for i1 = 1:numel(depths)
    for i2 = 1:numel(centers)
        for i3 = 1:numel(extents)
            r = center/2:0.01:(extent + center);
            f = @(x) potentialParametersFun(x,[depths(i1),centers(i2),extents(i3)],r);
            potentialParameters(i1,i2,i3,:) = fminsearch(f,[-depths(i1),centers(i2),(extents(i3)-centers(i2))/3],options);
        end
    end
    fprintf('%d/%d...\n',i1,numel(depths))
end
t = toc;
fprintf('time = %0.2f\naverage time = %0.5f\n',t, t/(numel(potentialParameters)/3))

%%
depth = -1;
center = 2;
extent = 10;

[~,i1] = min(abs(potentialParameters.depth-depth));
[~,i2] = min(abs(potentialParameters.center-center));
[~,i3] = min(abs(potentialParameters.extent-extent));

y = permute(potentialParameters.parameters(i1,i2,i3,:),[4,1,2,3]);

% r = center/2:0.01:(extent + center);
% 

% f = @(x) potentialParametersFun(x,[depth,center,extent],r);
% options = optimset('TolX',1e-5,'TolFun',1e-7,'MaxFunEvals',500);
% tic
% y = fminsearch(f,[-depth,center,(extent-center)/3],options);
% toc
% y
alpha = y(1);
cnt = y(2);
sig = y(3);


% %%


% f = @(x) min(Vint(r,mu,sig,x))-A);
% 
% options = optimset('Display','iter','TolX',1e-2);
% tic
% alpha = fzero(f,[0.1,2],options);
% toc
% alpha
r = 1e-2:0.01:20;
fig100 = figure(100);
clf(fig100);
V_exp = @(r,mu,sig,A) A*exp(-(r-mu).^2/(2*sig^2));
Vint = @(r,mu,sig,A) 1./r - A*exp(-(r-mu).^2/(2*sig^2));
dVint = @(D,mu,sig,A) -1./D.^2 + (A*(D-mu)/(sig^2)) .* exp(-(D-mu).^2/(2*sig^2));
G1 = V_exp(r,cnt,sig,alpha);

% tmp = 1./r - V_exp(r,2,2,0.5); [~,ind] = min(tmp); r(ind)

line(r,G1,'Color','g','LineWidth',2)
line(r,Vint(r,cnt,sig,alpha),'Color','y','LineWidth',2)
line(r,1./r,'Color','b','LineWidth',2)

title(sprintf('d=%0.2f, c=%0.1f, e=%0.0f',depth,center,extent))
legend({'gauss','1/r - gauss','1/r'})
ylim([-2,2])
drawnow;
line([0,20],[0 0],'Color',0.5*[1 1 1],'LineStyle','--')
goDark(gcf)
%%
figure
line(r(2:end-1),diff(1./r - G1,2))

%%
% 
% 
% close all
% 
% KAPPA_SMOOTHING_SIGMA = 2;
% 
% [Idsm,Idssm,BWsm] = getDeclumpTestCase(1);
% 
% % NOTE: speed up inpolygon inside getCurvatureAndShapeMarkers()
% [B,Bnorm,~,markers] = computeBoundaryInformation(BWsm,KAPPA_SMOOTHING_SIGMA); 
% 
% 
% B = B{1};
% Bnorm = Bnorm{1};
% markers = markers{1};
% %%
% [centers,Info] = computeObjectCenters(B,markers);
% 
% figure
% imshow(BWsm)
% hold on
% plot(centers(:,2),centers(:,1),'ro')
% 
% %%
% tmpFun = @() computeBoundaryInformation(BWsm,KAPPA_SMOOTHING_SIGMA);
% profile on
% timeit(tmpFun,4)
% profile off
% profile viewer
% 
% %%
% 
% BWtmp = createMaskFromBoundary(B);
% 
% tmpFun = @() createMaskFromBoundary(B);
% 
% timeit(tmpFun,1)
% 
% 
% %%
% figure
% imshow(BWtmp)
% 
% %%
% figure
% imshow(BWsm)
% hold on
% plot(markers(:,2),markers(:,1),'b.')
% quiver(B(:,2),B(:,1),Bnorm(:,2),Bnorm(:,1),0.25,'color','r')