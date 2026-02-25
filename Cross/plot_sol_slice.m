clc,clear
close all

xa=0;xb=1;yc=0;yd=1;

%导入参考解
Data_r = load('Data_r.mat');

%依次导入3^i,i=2,...,5解
Data_3e2 = load('Data_3e2.mat');
Data_3e3 = load('Data_3e3.mat');
Data_3e4 = load('Data_3e4.mat');
Data_3e5 = load('Data_3e5.mat');

%绘制3^2剖分下数值解图和3^5剖分下数值解图
%3^2
nx1 = sqrt(Data_3e2.dof);
ny1 = nx1;
[x_rect1, y_rect1] = meshgrid(linspace(xa, xb, nx1), linspace(yc, yd, ny1));
node1=[x_rect1(:),y_rect1(:)];
X_meshgrid = reshape(node1(:,1),nx1,ny1);
Y_meshgrid = reshape(node1(:,2),nx1,ny1);
U_meshgrid = reshape(Data_3e2.uc,nx1,ny1);
Interp_N = nx1; Interp_M = ny1; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq1 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq1)
colormap('parula');xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_3e2_Solution','epsc');

% %3^3
% nx2 = sqrt(Data_3e3.dof);
% ny2 = nx2;
% [x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
% node2=[x_rect2(:),y_rect2(:)];
% X_meshgrid = reshape(node2(:,1),nx2,ny2);
% Y_meshgrid = reshape(node2(:,2),nx2,ny2);
% U_meshgrid = reshape(Data_3e3.uc,nx2,ny2);
% Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
% [Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
% Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
% figure; 
% hold on;
% surf(Xq,Yq,Vq2)
% colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
% contour(X_meshgrid,Y_meshgrid,U_meshgrid')
% hold off;
% saveas(gcf, 'Cross_3e3_Solution','epsc');
% 
% %3^4
% nx2 = sqrt(Data_3e4.dof);
% ny2 = nx2;
% [x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
% node2=[x_rect2(:),y_rect2(:)];
% X_meshgrid = reshape(node2(:,1),nx2,ny2);
% Y_meshgrid = reshape(node2(:,2),nx2,ny2);
% U_meshgrid = reshape(Data_3e4.uc,nx2,ny2);
% Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
% [Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
% Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
% figure; 
% hold on;
% surf(Xq,Yq,Vq2)
% colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
% contour(X_meshgrid,Y_meshgrid,U_meshgrid')
% hold off;
% saveas(gcf, 'Cross_3e4_Solution','epsc');
% 
% %3^5
% nx2 = sqrt(Data_3e5.dof);
% ny2 = nx2;
% [x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
% node2=[x_rect2(:),y_rect2(:)];
% X_meshgrid = reshape(node2(:,1),nx2,ny2);
% Y_meshgrid = reshape(node2(:,2),nx2,ny2);
% U_meshgrid = reshape(Data_3e5.uc,nx2,ny2);
% Interp_N = 121; Interp_M = 121; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
% [Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
% Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
% figure; 
% hold on;
% surf(Xq,Yq,Vq2)
% colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
% contour(X_meshgrid,Y_meshgrid,U_meshgrid')
% hold off;
% saveas(gcf, 'Cross_3e5_Solution','epsc');
% 
% %slice
% %y=0.45
% figure;
% plot(Data_r.R_slice0_51,Data_r.U_slice0_51,'LineWidth', 1.5);
% hold on;
% plot(Data_3e2.R_slice0_51,Data_3e2.U_slice0_51,'LineWidth', 1.5);
% plot(Data_3e3.R_slice0_51,Data_3e3.U_slice0_51,'LineWidth', 1.5);
% plot(Data_3e4.R_slice0_51,Data_3e4.U_slice0_51,'LineWidth', 1.5);
% plot(Data_3e5.R_slice0_51,Data_3e5.U_slice0_51,'LineWidth', 1.5);
% hold off;
% legend('Reference','3^2 \times 3^2 mesh','3^3 \times 3^3 mesh','3^4 \times 3^4 mesh', '3^5 \times 3^5 mesh', 'Location', 'northeast'); 
% saveas(gcf, 'Cross_Slice_y_0p45','epsc'); 
% 
% %y=0.5
% figure;
% plot(Data_r.R_slice1_51,Data_r.U_slice1_51,'LineWidth', 1.5);
% hold on;
% plot(Data_3e2.R_slice1_51,Data_3e2.U_slice1_51,'LineWidth', 1.5);
% plot(Data_3e3.R_slice1_51,Data_3e3.U_slice1_51,'LineWidth', 1.5);
% plot(Data_3e4.R_slice1_51,Data_3e4.U_slice1_51,'LineWidth', 1.5);
% plot(Data_3e5.R_slice1_51,Data_3e5.U_slice1_51,'LineWidth', 1.5);
% hold off;
% legend('Reference','3^2 \times 3^2 mesh','3^3 \times 3^3 mesh','3^4 \times 3^4 mesh', '3^5 \times 3^5 mesh', 'Location', 'northeast'); 
% saveas(gcf, 'Cross_Slice_y_0p5','epsc'); 
% 
% %x=0.3
% figure;
% plot(Data_r.R_slice2_51,Data_r.U_slice2_51,'LineWidth', 1.5);
% hold on;
% plot(Data_3e2.R_slice2_51,Data_3e2.U_slice2_51,'LineWidth', 1.5);
% plot(Data_3e3.R_slice2_51,Data_3e3.U_slice2_51,'LineWidth', 1.5);
% plot(Data_3e4.R_slice2_51,Data_3e4.U_slice2_51,'LineWidth', 1.5);
% plot(Data_3e5.R_slice2_51,Data_3e5.U_slice2_51,'LineWidth', 1.5);
% hold off;
% legend('Reference','3^2 \times 3^2 mesh','3^3 \times 3^3 mesh','3^4 \times 3^4 mesh', '3^5 \times 3^5 mesh', 'Location', 'northeast');
% saveas(gcf, 'Cross_Slice_x_0p3','epsc'); 
% 
% %x=0.4
% figure;
% plot(Data_r.R_slice3_51,Data_r.U_slice3_51,'LineWidth', 1.5);
% hold on;
% plot(Data_3e2.R_slice3_51,Data_3e2.U_slice3_51,'LineWidth', 1.5);
% plot(Data_3e3.R_slice3_51,Data_3e3.U_slice3_51,'LineWidth', 1.5);
% plot(Data_3e4.R_slice3_51,Data_3e4.U_slice3_51,'LineWidth', 1.5);
% plot(Data_3e5.R_slice3_51,Data_3e5.U_slice3_51,'LineWidth', 1.5);
% hold off;
% legend('Reference','3^2 \times 3^2 mesh','3^3 \times 3^3 mesh','3^4 \times 3^4 mesh', '3^5 \times 3^5 mesh', 'Location', 'northeast');
% saveas(gcf, 'Cross_Slice_x_0p4','epsc');
% % 对比图
% 
% 
% %% 以下为废弃调色
% 
% % 'Color', [0.2 0.4 0.8 1.0];
% % 'Color', [1.0 0.5 0.1 0.7];
% % 'Color', [0.6 0.2 0.7 0.7];
% % 'Color', [0.2 0.7 0.3 0.7];
% % 'Color', [0.8 0.2 0.2 0.7];
% 
