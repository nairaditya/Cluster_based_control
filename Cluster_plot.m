function Cluster_plot(nclust,CL,CD,CL_dot,c1_Centers,c1_Labels,ai)

%% rescale to original values

minx_a     = min(CL);
denx_a     = max(CL) - min(CL);
miny_a     = min(CL_dot);
deny_a     = max(CL_dot) - min(CL_dot);
minz_a     = min(CD);
denz_a     = max(CD) - min(CD);
ai(:,1)    = ai(:,1)*denx_a+minx_a;
ai(:,2)    = ai(:,2)*deny_a+miny_a;
ai(:,3)    = ai(:,3)*denz_a+minz_a;
c1_Centers(:,1) = c1_Centers(:,1)*denx_a+minx_a;
c1_Centers(:,2) = c1_Centers(:,2)*deny_a+miny_a;
c1_Centers(:,3) = c1_Centers(:,3)*denz_a+minz_a;

%% Parameters for plotting clusters

nCluster1 = length(unique(c1_Labels));
nCluster  = nCluster1;
CMap      = jet(nCluster);
step      = floor(nCluster/nCluster1);
CMap      = CMap(1:step:end,:);
LineWidth = 2;
TextSize  = 12;
MarkerSize= 8;

figure;subplot(111);
colormap(CMap);

%% Determine Axes

xMin = min(min([ai(:,1)]));
xMax = max(max([ai(:,1)]));
yMin = min(min([ai(:,2)]));
yMax = max(max([ai(:,2)]));
zMin = min(min([ai(:,3)]));
zMax = max(max([ai(:,3)]));
dx   = 0.05*(xMax-xMin);
dy   = 0.05*(yMax-yMin);
dz   = 0.05*(zMax-zMin);
xMin = xMin - dx;
xMax = xMax + dx;
yMin = yMin - dy;
yMax = yMax + dy;
zMin = zMin - dz;
zMax = zMax + dz;

%% Plot data points 

hold on;
for iCluster = 1:nCluster1 
    plot3(ai(c1_Labels==iCluster,1),ai(c1_Labels==iCluster,2),...
        ai(c1_Labels==iCluster,3),'.','MarkerSize',2*MarkerSize,...
        'MarkerEdgeColor',CMap(iCluster,:),'Color',CMap(iCluster,:),...
        'LineWidth',LineWidth)
end

%% Plot centroids

for iCluster = 1:size(c1_Centers,1)
    h(iCluster) = plot3(c1_Centers(iCluster,1), c1_Centers(iCluster,2), ...
        c1_Centers(iCluster,3),'s','MarkerEdgeColor','k',...
        'MarkerFaceColor',CMap(iCluster,:),'MarkerSize',2*MarkerSize,...
        'LineWidth',LineWidth);
end

centroids_txt = [1:size(c1_Centers,1)];
for iCluster = 1:size(c1_Centers,1)
    text(c1_Centers(iCluster,1)+dx*100*abs(c1_Centers(iCluster,1)), ...
        c1_Centers(iCluster,2)+dy*100*abs(c1_Centers(iCluster,2)),...
        c1_Centers(iCluster,3)+dz*100*abs(c1_Centers(iCluster,3)),...
        [num2str(centroids_txt(iCluster))],'FontSize',TextSize);
end
view(3);
grid on;box off;
ylabel('$\dot{C_L}$','Interpreter','latex','Fontsize',TextSize+2);
xlabel('$C_L$','Interpreter','latex','Fontsize',TextSize+2);
zlabel('$C_D$','Interpreter','latex','Fontsize',TextSize+2);
axis tight;
set(gca, 'Fontsize', TextSize+2);
print('-dpng',['clusterCL',num2str(nclust),'.png']);
close;

