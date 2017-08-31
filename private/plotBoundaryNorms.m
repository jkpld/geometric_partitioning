function plotBoundaryNorms(B,N)

figure
hold on
plot(B(:,2), B(:,1),'Color','k','LineWidth',2)
quiver(B(:,2),B(:,1),N(:,2),N(:,1),0.25,'Color',0.5*[1 1 1],'LineWidth',1.5);

end