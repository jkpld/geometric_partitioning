function [fig,ax] = plotObjectCutDebugInfo(B,N,cuts,Info)
% PLOTOBJECTCUTDEBUGINFO  Helper function for plotting debuggin information

% James Kapaldo
% 2016-10-29
Info.topLeft = [0,0];
r0 = [nan, nan];
r_end = [nan, nan];
centers = Info.centers - Info.topLeft;
% r0 = Info.r0 - Info.topLeft;
% r_end = Info.r_end - Info.topLeft;

fig = figure;
ax = axes('Parent',fig);
hold on

aC = Info.CreateObjectPartition.assignedCenter;
assignmentLine = nan(3*size(B,1),2);
ind = 1:3:3*size(B,1);
assignmentLine(ind,:) = B;
assignmentLine(ind(~isnan(aC))+1,:) = centers(aC(~isnan(aC)),:);

p1 = plot(assignmentLine(:,1), assignmentLine(:,2),'Color',0.4*[1 1 1]);
plot(B(:,1), B(:,2),'Color','k','LineWidth',2)
p2 = plot(centers(:,1),centers(:,2),'bo','MarkerFaceColor','b');
p3 = quiver(B(:,1),B(:,2),N(:,1),N(:,2),0.25,'Color',0.5*[1 1 1],'LineWidth',1.5);
p4 = plot(r0(:,1),r0(:,2),'r.');
p5 = plot(r_end(:,1),r_end(:,2),'b.');

for i = 1:numel(Info.CreateObjectPartition.CalculateCuts.partitions)
    pP = line(Info.CreateObjectPartition.CalculateCuts.partitions{i}(:,1),Info.CreateObjectPartition.CalculateCuts.partitions{i}(:,2),'color','g');
end

if isempty(cuts)
    cuts = [nan,nan,nan,nan];
end
for i = 1:size(cuts,1)
    tmp = (cuts(i,1:2) + cuts(i,3:4))/2;
    pC = line(cuts(i,[1,3]),cuts(i,[2,4]),'color',[0 0 1],'lineWidth',2);
    text(tmp(1),tmp(2),num2str(i),'FontSize',15,'fontweight','bold','color',[0 0 0.5],'Tag','ignore','Clipping','on')
end

p6 = plot(cuts(:,1),cuts(:,2),'o','LineWidth',1,'Color',0.2*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g');
p7 = plot(cuts(:,3),cuts(:,4),'d','LineWidth',1,'Color',0.2*[1 1 1],'MarkerSize',7,'MarkerFaceColor','g');

legend([p1,p3,p4,p5,p2,pP,pC,p6,p7],{'Boundary assignments','Boundary normals','Initial particle locations','Final particle locations','Centers','Partitions before optimization','Cuts','Cut starts','Cut ends'},'Location','eastoutside')

for i = 1:size(centers,1)
    text(centers(i,1),centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color',[0.5 0 0],'Tag','ignore','Clipping','on')
end

grid on
box on
axis tight
daspect([1 1 1])
drawnow;


end
