%% compare 
clc,clear
close all

%两个参考解
Ref_1 = load('referenceData.mat');
Ref_2 = load('referenceData2.mat');

%圆形裂缝解
curve_20_sol1 = load('Curve_20_1.mat');
curve_40_sol1 = load('Curve_40_1.mat');
curve_80_sol1 = load('Curve_80_1.mat');

%半圆形裂缝解
curve_20_sol2 = load('Curve_20_2.mat');
curve_40_sol2 = load('Curve_40_2.mat');
curve_80_sol2 = load('Curve_80_2.mat');

%圆形裂缝对比参考解
figure; % slice 1 y=0.5
hold on;
plot(Ref_1.R_slice1,Ref_1.U_slice1,'LineWidth',1);
plot(curve_20_sol1.R_slice1,curve_20_sol1.U_slice1,'LineWidth',1);
plot(curve_40_sol1.R_slice1,curve_40_sol1.U_slice1,'LineWidth',1);
plot(curve_80_sol1.R_slice1,curve_80_sol1.U_slice1,'LineWidth',1);
legend('Reference','20 \times 20 mesh','40 \times 40 mesh','80 \times 80 mesh');
hold off;
saveas(gcf, 'Curve_Slice_y_0p5_1','epsc'); 


figure; % slice 1 x=0.4
hold on;
plot(Ref_1.R_slice3,Ref_1.U_slice3,'LineWidth',1);
plot(curve_20_sol1.R_slice2,curve_20_sol1.U_slice2,'LineWidth',1);
plot(curve_40_sol1.R_slice2,curve_40_sol1.U_slice2,'LineWidth',1);
plot(curve_80_sol1.R_slice2,curve_80_sol1.U_slice2,'LineWidth',1);
legend('Reference','20 \times 20 mesh','40 \times 40 mesh','80 \times 80 mesh');
hold off;
saveas(gcf, 'Curve_Slice_x_0p4_1','epsc'); 

%半圆形裂缝对比参考解
figure; % slice 1 y=0.5
hold on;
plot(Ref_2.R_slice1,Ref_2.U_slice1,'LineWidth',1);
plot(curve_20_sol2.R_slice1,curve_20_sol2.U_slice1,'LineWidth',1);
plot(curve_40_sol2.R_slice1,curve_40_sol2.U_slice1,'LineWidth',1);
plot(curve_80_sol2.R_slice1,curve_80_sol2.U_slice1,'LineWidth',1);
legend('Reference','20 \times 20 mesh','40 \times 40 mesh','80 \times 80 mesh');
hold off;
saveas(gcf, 'Curve_Slice_y_0p5_2','epsc'); 


figure; % slice 1 x=0.4
hold on;
plot(Ref_2.R_slice3,Ref_2.U_slice3,'LineWidth',1);
plot(curve_20_sol2.R_slice2,curve_20_sol2.U_slice2,'LineWidth',1);
plot(curve_40_sol2.R_slice2,curve_40_sol2.U_slice2,'LineWidth',1);
plot(curve_80_sol2.R_slice2,curve_80_sol2.U_slice2,'LineWidth',1);
legend('Reference','20 \times 20 mesh','40 \times 40 mesh','80 \times 80 mesh');
hold off;
saveas(gcf, 'Curve_Slice_x_0p4_2','epsc'); 
