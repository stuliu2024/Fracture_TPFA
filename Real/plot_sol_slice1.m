clc,clear
close all


xa=0;xb=700;yc=0;yd=600;

addpath(genpath(pwd))

Data_105 = load('Data_560.mat');

% compare1
Data_ldg1 = load('Data_ldg_350.mat');

% ldg
Data_ldg1.Slice_Res = mat2cell(Data_ldg1.Slice_Res, Data_ldg1.length_slice, 1);
Data_ldg1.Slice_Ues = mat2cell(Data_ldg1.Slice_Ues, Data_ldg1.length_slice, 1);

Slice_Res_ldg11 = Data_ldg1.Slice_Res{1, 1};
Slice_Res_ldg21 = Data_ldg1.Slice_Res{2, 1};

Slice_Ues_ldg11 = Data_ldg1.Slice_Ues{1, 1};
Slice_Ues_ldg21 = Data_ldg1.Slice_Ues{2, 1};


%cg
Data_cg_105 = load('Data.mat');
% Data_cg_175 = load('Data_cg_175.mat');

%slice
%y=500
figure;
plot(Data_105.R_slice1_51,Data_105.U_slice1_51,'LineWidth', 1.5);
hold on;
plot(Data_cg_105.R_slice1_51,Data_cg_105.U_slice1_51,'LineWidth', 1.5);
plot(Slice_Res_ldg11,Slice_Ues_ldg11,'LineWidth', 1.5);
hold off;
legend('RDFM-TPFA','RDFM-CG','RDFM-LDG', 'Location', 'northeast');   
saveas(gcf, 'Real_Slice_y_500_compare1','epsc'); 

%x=625
figure;
plot(Data_105.R_slice1_52,Data_105.U_slice1_52,'LineWidth', 1.5);
hold on;
plot(Data_cg_105.R_slice2_51,Data_cg_105.U_slice2_51,'LineWidth', 1.5);
plot(Slice_Res_ldg21,Slice_Ues_ldg21,'LineWidth', 1.5);
hold off;
legend('RDFM-TPFA','RDFM-CG','RDFM-LDG', 'Location', 'northeast');     
saveas(gcf, 'Real_Slice_x_625_compare1','epsc'); 
