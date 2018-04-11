function plot_partition_results(BWpart, r0, region)

R = region;
if iscell(r0)
    r0 = cat(1,r0{:});
end

L = bwlabel(BWpart(R(1,1):R(2,1),R(1,2):R(2,2)));

cmap = [0.87843      0.95294      0.85882;
      0.65882      0.86667       0.7098;
      0.30588      0.70196      0.82745;
     0.031373      0.40784      0.67451];
 
cmap = repmat(cmap,[ceil(max(L(:))/4),1]);


figure
imshow(label2rgb(L,cmap,'k','shuffle'))
hold on;
valid = r0(:,2)>R(1,2) & r0(:,2)<R(2,2) & r0(:,1)>R(1,1) & r0(:,1)<R(1,2);

plot(r0(valid,2)-R(1,2),r0(valid,1)-R(1,1),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',3,'LineWidth',0.5)

title('Partitioning results')
ax = gca;
ax.XRuler.Visible = 'off';
ax.YRuler.Visible = 'off';
setTheme(gcf,'light')
end