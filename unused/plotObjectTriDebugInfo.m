function [fig,ax] = plotObjectTriDebugInfo(B,cut,cutK,Info)
% PLOTOBJECTTRIDEBUGINFO  Helper function for plotting debuggin information
% about triangle cuts

% James Kapaldo
% 2016-10-29

centers = Info.centers - Info.topLeft;
TriInfo = Info.CreateTriangleCuts;

fig = figure;
ax = axes('Parent',fig);
hold on

plot(B(:,1), B(:,2),'Color','k','LineWidth',2)
pC = plot(centers(:,1),centers(:,2),'bo','MarkerFaceColor','b');

for group = 1:numel(TriInfo)
    
    for i = 1:size(TriInfo(group).triGroup,2)
        triCent = TriInfo(group).triGroup(:,i);
        p4 = plot(centers(triCent([1,2]),1),centers(triCent([1,2]),2),'Color',0.3*[1 1 1],'LineWidth',2);
        plot(centers(triCent([1,3]),1),centers(triCent([1,3]),2),'Color',0.3*[1 1 1],'LineWidth',2);
        plot(centers(triCent([2,3]),1),centers(triCent([2,3]),2),'Color',0.3*[1 1 1],'LineWidth',2);
    end
    
    pM = plot(TriInfo(group).triangleCenterSearchPoints(:,1)-Info.topLeft(:,1),TriInfo(group).triangleCenterSearchPoints(:,2)-Info.topLeft(:,2),'ys','MarkerSize',4,'LineWidth',0.5);
    for i = 1:size(TriInfo(group).originalCuts)
        p1 = plot(TriInfo(group).originalCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).originalCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[1 0 0],'LineWidth',2);
        p2 = plot(TriInfo(group).vertexOptimizedCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).vertexOptimizedCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[0 1 0],'LineWidth',2);
        p3 = plot(TriInfo(group).optimizedCuts(i,[1,3])-Info.topLeft(:,1), TriInfo(group).optimizedCuts(i,[2,4])-Info.topLeft(:,2),'Color',0.5*[0 0 1],'LineWidth',2);
    end
    
    cutG = cut{group}(isnan(cutK{group}(:,1)),:);
    for i = 1:size(cutG,1)
        plot(cutG(i,[1,3]),cutG(i,[2,4]),'Color',0.5*[0 0 1],'LineWidth',2,'LineStyle','--')
    end
end

legend([pC,p4,pM,p1,p2,p3],{'Centers','Triangle edges','Triangle center search points','Original cuts','Vertex optimized cuts','OptimizedCuts'},'Location','eastoutside')

for i = 1:size(centers,1)
    text(centers(i,1),centers(i,2),num2str(i),'FontSize',15,'fontweight','bold','color',[0.5 0 0],'Tag','ignore','Clipping','on')
end

grid on
box on
axis tight
daspect([1 1 1])
drawnow;

end