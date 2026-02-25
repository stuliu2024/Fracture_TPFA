clc,clear
close all

xa=0;xb=1;yc=0;yd=1;

%导入参考解
Data_r = load('Data_r.mat');

%依次导入3^i,i=2,...,5解
Data_20 = load('Data_20.mat');
Data_30 = load('Data_30.mat');
Data_60 = load('Data_60.mat');

%绘制3^2剖分下数值解图和3^5剖分下数值解图
%10
nx1 = sqrt(Data_20.dof);
ny1 = nx1;
[x_rect1, y_rect1] = meshgrid(linspace(xa, xb, nx1), linspace(yc, yd, ny1));
node1=[x_rect1(:),y_rect1(:)];
X_meshgrid = reshape(node1(:,1),nx1,ny1);
Y_meshgrid = reshape(node1(:,2),nx1,ny1);
U_meshgrid = reshape(Data_20.uc,nx1,ny1);
Interp_N = nx1; Interp_M = ny1; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq1 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq1)
colormap('parula');xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_20_Solution','epsc');

%40
nx2 = sqrt(Data_30.dof);
ny2 = nx2;
[x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
node2=[x_rect2(:),y_rect2(:)];
X_meshgrid = reshape(node2(:,1),nx2,ny2);
Y_meshgrid = reshape(node2(:,2),nx2,ny2);
U_meshgrid = reshape(Data_30.uc,nx2,ny2);
Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq2)
colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_30_Solution','epsc');

%80
nx2 = sqrt(Data_60.dof);
ny2 = nx2;
[x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
node2=[x_rect2(:),y_rect2(:)];
X_meshgrid = reshape(node2(:,1),nx2,ny2);
Y_meshgrid = reshape(node2(:,2),nx2,ny2);
U_meshgrid = reshape(Data_60.uc,nx2,ny2);
Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq2)
colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_60_Solution','epsc');


%slice
%y=0.45
figure;
plot(Data_r.R_slice0_51,Data_r.U_slice0_51,'LineWidth', 1.5);
hold on;
plot(Data_20.R_slice0_51,Data_20.U_slice0_51,'LineWidth', 1.5);
plot(Data_30.R_slice0_51,Data_30.U_slice0_51,'LineWidth', 1.5);
plot(Data_60.R_slice0_51,Data_60.U_slice0_51,'LineWidth', 1.5);
hold off;
legend('Reference','20 \times 20 mesh','30 \times 30 mesh','60 \times 60 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_y_0p45_1','epsc'); 



%x=0.3
figure;
plot(Data_r.R_slice2_51,Data_r.U_slice2_51,'LineWidth', 1.5);
hold on;
plot(Data_20.R_slice2_51,Data_20.U_slice2_51,'LineWidth', 1.5);
plot(Data_30.R_slice2_51,Data_30.U_slice2_51,'LineWidth', 1.5);
plot(Data_60.R_slice2_51,Data_60.U_slice2_51,'LineWidth', 1.5);
hold off;
legend('Reference','20 \times 20 mesh','30 \times 30 mesh','60 \times 60 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_x_0p3_1','epsc'); 

%x=0.4
figure;
plot(Data_r.R_slice3_51,Data_r.U_slice3_51,'LineWidth', 1.5);
hold on;
plot(Data_20.R_slice3_51,Data_20.U_slice3_51,'LineWidth', 1.5);
plot(Data_30.R_slice3_51,Data_30.U_slice3_51,'LineWidth', 1.5);
plot(Data_60.R_slice3_51,Data_60.U_slice3_51,'LineWidth', 1.5);
hold off;
legend('Reference','20 \times 20 mesh','30 \times 30 mesh','60 \times 60 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_x_0p4_1','epsc');



%%-----------------------------------odd---------------------------%
%依次导入3^i,i=2,...,5解
Data_21 = load('Data_21.mat');
Data_31 = load('Data_31.mat');
Data_61 = load('Data_61.mat');

%绘制3^2剖分下数值解图和3^5剖分下数值解图
%10
nx1 = sqrt(Data_21.dof);
ny1 = nx1;
[x_rect1, y_rect1] = meshgrid(linspace(xa, xb, nx1), linspace(yc, yd, ny1));
node1=[x_rect1(:),y_rect1(:)];
X_meshgrid = reshape(node1(:,1),nx1,ny1);
Y_meshgrid = reshape(node1(:,2),nx1,ny1);
U_meshgrid = reshape(Data_21.uc,nx1,ny1);
Interp_N = nx1; Interp_M = ny1; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq1 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq1)
colormap('parula');xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_21_Solution','epsc');

%40
nx2 = sqrt(Data_31.dof);
ny2 = nx2;
[x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
node2=[x_rect2(:),y_rect2(:)];
X_meshgrid = reshape(node2(:,1),nx2,ny2);
Y_meshgrid = reshape(node2(:,2),nx2,ny2);
U_meshgrid = reshape(Data_31.uc,nx2,ny2);
Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq2)
colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_31_Solution','epsc');

%80
nx2 = sqrt(Data_61.dof);
ny2 = nx2;
[x_rect2, y_rect2] = meshgrid(linspace(xa, xb, nx2), linspace(yc, yd, ny2));
node2=[x_rect2(:),y_rect2(:)];
X_meshgrid = reshape(node2(:,1),nx2,ny2);
Y_meshgrid = reshape(node2(:,2),nx2,ny2);
U_meshgrid = reshape(Data_61.uc,nx2,ny2);
Interp_N = nx2; Interp_M = ny2; Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
[Xq,Yq] = meshgrid(linspace(Interp_xmin,Interp_xmax,Interp_N),linspace(Interp_ymin,Interp_ymax,Interp_M)); 
Vq2 = interp2(X_meshgrid,Y_meshgrid,U_meshgrid',Xq,Yq,'linear');
figure; 
hold on;
surf(Xq,Yq,Vq2)
colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);view(30,30);
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
hold off;
saveas(gcf, 'Cross_61_Solution','epsc');


%slice
%y=0.45
figure;
plot(Data_r.R_slice0_51,Data_r.U_slice0_51,'LineWidth', 1.5);
hold on;
plot(Data_21.R_slice0_51,Data_21.U_slice0_51,'LineWidth', 1.5);
plot(Data_31.R_slice0_51,Data_31.U_slice0_51,'LineWidth', 1.5);
plot(Data_61.R_slice0_51,Data_61.U_slice0_51,'LineWidth', 1.5);
hold off;
legend('Reference','21 \times 21 mesh','31 \times 31 mesh','61 \times 61 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_y_0p45_2','epsc'); 



%x=0.3
figure;
plot(Data_r.R_slice2_51,Data_r.U_slice2_51,'LineWidth', 1.5);
hold on;
plot(Data_21.R_slice2_51,Data_21.U_slice2_51,'LineWidth', 1.5);
plot(Data_31.R_slice2_51,Data_31.U_slice2_51,'LineWidth', 1.5);
plot(Data_61.R_slice2_51,Data_61.U_slice2_51,'LineWidth', 1.5);
hold off;
legend('Reference','21 \times 21 mesh','31 \times 31 mesh','61 \times 61 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_x_0p3_2','epsc'); 

%x=0.4
figure;
plot(Data_r.R_slice3_51,Data_r.U_slice3_51,'LineWidth', 1.5);
hold on;
plot(Data_21.R_slice3_51,Data_21.U_slice3_51,'LineWidth', 1.5);
plot(Data_31.R_slice3_51,Data_31.U_slice3_51,'LineWidth', 1.5);
plot(Data_61.R_slice3_51,Data_61.U_slice3_51,'LineWidth', 1.5);
hold off;
legend('Reference','21 \times 21 mesh','31 \times 31 mesh','61 \times 61 mesh', 'Location', 'northeast'); 
saveas(gcf, 'Cross_Slice_x_0p4_2','epsc');


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
