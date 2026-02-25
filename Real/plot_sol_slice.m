clc,clear
close all

xa=0;xb=700;yc=0;yd=600;

%导入参考解
% Data_r = load('mfd_geiger_conductive.mat');

%依次导入3^i,i=2,...,5解
Data_r = load('Data.mat');
Data_105 = load('Data_105.mat');
Data_175 = load('Data_175.mat');
Data_560 = load('Data_560.mat');
%slice
%y=500
figure;
plot(Data_r.R_slice1_51,Data_r.U_slice1_51,'LineWidth', 1.5);
hold on;
plot(Data_105.R_slice1_51,Data_105.U_slice1_51,'LineWidth', 1.5);
plot(Data_175.R_slice1_51,Data_175.U_slice1_51,'LineWidth', 1.5);
hold off;
legend('Reference','105 \times 90 mesh','175 \times 150 mesh', 'Location', 'northeast');   
saveas(gcf, 'Real_Slice_y_500','epsc'); 

%x=625
figure;
plot(Data_r.R_slice2_51,Data_r.U_slice2_51,'LineWidth', 1.5);
hold on;
plot(Data_105.R_slice1_52,Data_105.U_slice1_52,'LineWidth', 1.5);
plot(Data_175.R_slice1_52,Data_175.U_slice1_52,'LineWidth', 1.5);
hold off;
legend('Reference','105 \times 90 mesh','175 \times 150 mesh', 'Location', 'northeast');   
saveas(gcf, 'Real_Slice_x_625','epsc'); 



% compare
Data_ldg = load('Data_ldg_175.mat');

% ldg
Data_ldg.Slice_Res = mat2cell(Data_ldg.Slice_Res, Data_ldg.length_slice, 1);
Data_ldg.Slice_Ues = mat2cell(Data_ldg.Slice_Ues, Data_ldg.length_slice, 1);

Slice_Res_ldg1 = Data_ldg.Slice_Res{1, 1};
Slice_Res_ldg2 = Data_ldg.Slice_Res{2, 1};

Slice_Ues_ldg1 = Data_ldg.Slice_Ues{1, 1};
Slice_Ues_ldg2 = Data_ldg.Slice_Ues{2, 1};



%slice
%y=500
figure;
plot(Data_r.R_slice1_51,Data_r.U_slice1_51,'LineWidth', 1.5);
hold on;
plot(Data_560.R_slice1_51,Data_560.U_slice1_51,'LineWidth', 1.5);
plot(Slice_Res_ldg1,Slice_Ues_ldg1,'LineWidth', 1.5);
hold off;
legend('Reference','RDFM-TPFA','RDFM-LDG', 'Location', 'northeast');   
saveas(gcf, 'Real_Slice_y_500_compare','epsc'); 

%x=625
figure;
plot(Data_r.R_slice2_51,Data_r.U_slice2_51,'LineWidth', 1.5);
hold on;
plot(Data_560.R_slice1_52,Data_560.U_slice1_52,'LineWidth', 1.5);
plot(Slice_Res_ldg2,Slice_Ues_ldg2,'LineWidth', 1.5);
hold off;
legend('Reference','RDFM-TPFA','RDFM-LDG', 'Location', 'northeast');   
saveas(gcf, 'Real_Slice_x_625_compare','epsc'); 
