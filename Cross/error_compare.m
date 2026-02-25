clc;
clear;
close all;

Data_r = load('Data_r.mat');
Data_485 = load('Data_121.mat');
Data_cg = load('Data_cg.mat');
Data_ldg = load('Data_ldg.mat');

%Cpu_time
Cpu_time_tpfa = Data_485.Cpu_time
Cpu_time_cg = Data_cg.Cpu_time
Cpu_time_ldg = Data_ldg.Cpu_time

%Dof
dof_tpfa = Data_485.dof
dof_cg = Data_cg.dof
dof_ldg = Data_ldg.dof

%condition numbers
conditions_tpfa = Data_485.conditionalnumber
conditions_cg = Data_cg.conditionalnumber

%sparsity
sparsity_tpfa = Data_485.sparsity
sparsity_cg = Data_cg.sparsity


% %Matrix Errors && Fracture Errors
% %Tpfa
[totalRelativeMatrixError_tpfa] = compute_2d_error( 'Data_485.mat', 'Data_rrr.mat')
[totalRelativeFractureError_tpfa] = compute_1d_error( 'Data_485.mat', 'Data_rrr.mat')
% 
% %cg
% [totalRelativeMatrixError_cg] = compute_2d_error( 'Data_cg.mat', 'Data_rrr.mat')
% [totalRelativeFractureError_cg] = compute_1d_error( 'Data_cg.mat', 'Data_rrr.mat')
% 
% %ldg
% [totalRelativeMatrixError_ldg] = compute_2d_error( 'Data_ldg.mat', 'Data_rrr.mat')
% [totalRelativeFractureError_ldg] = compute_1d_error( 'Data_ldg.mat', 'Data_rrr.mat')

%slice

%ldg

% Data_ldg.Slice_Res = mat2cell(Data_ldg.Slice_Res, Data_ldg.length_slice, 1);
% Data_ldg.Slice_Ues = mat2cell(Data_ldg.Slice_Ues, Data_ldg.length_slice, 1);
% 
% Slice_Res_ldg1 = Data_ldg.Slice_Res{1, 1};
% Slice_Res_ldg2 = Data_ldg.Slice_Res{2, 1};
% Slice_Res_ldg3 = Data_ldg.Slice_Res{3, 1};
% Slice_Res_ldg4 = Data_ldg.Slice_Res{4, 1};
% 
% Slice_Ues_ldg1 = Data_ldg.Slice_Ues{1, 1};
% Slice_Ues_ldg2 = Data_ldg.Slice_Ues{2, 1};
% Slice_Ues_ldg3 = Data_ldg.Slice_Ues{3, 1};
% Slice_Ues_ldg4 = Data_ldg.Slice_Ues{4, 1};
% 
% % y=0.45
% figure;
% plot(Data_r.R_slice0_51, Data_r.U_slice0_51, 'LineWidth', 1.5);
% hold on;
% plot(Data_485.R_slice0_51, Data_485.U_slice0_51, 'LineWidth', 1.5);
% plot(Data_cg.R_slice0_51, Data_cg.U_slice0_51, 'LineWidth', 1.5);
% plot(Slice_Res_ldg1, Slice_Ues_ldg1, 'LineWidth', 1.5);
% hold off;
% legend('Reference', 'RDFM-TPFA', 'RDFM-CG', 'RDFM-LDG', 'Location', 'northeast');
% %局部图
% zoom_x_start = 0.5858;  % x1
% zoom_x_end = 0.5872;    % x2
% zoom_y_start = 0.4800;  % y1
% zoom_y_end = 0.4825;    % y2
% 
% % hold on;
% % rectangle('Position', [zoom_x_start, zoom_y_start, zoom_x_end - zoom_x_start, zoom_y_end - zoom_y_start],'EdgeColor', 'r', 'LineWidth', 3.5, 'LineStyle', '-');
% % hold off;
% inset_ax = axes('Position', [0.15, 0.15, 0.3, 0.3]);
% box on;
% grid off;
% set(gca, 'XTick',[],'YTick',[]);
% hold on;
% set(inset_ax, 'FontSize', 8);
% plot(inset_ax, Data_r.R_slice0_51, Data_r.U_slice0_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_485.R_slice0_51, Data_485.U_slice0_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_cg.R_slice0_51, Data_cg.U_slice0_51, 'LineWidth', 1.5);
% plot(inset_ax, Slice_Res_ldg1, Slice_Ues_ldg1, 'LineWidth', 1.5);
% hold(inset_ax, 'off');
% xlim(inset_ax, [zoom_x_start, zoom_x_end]);
% ylim(inset_ax, [zoom_y_start, zoom_y_end]);
% % annotation('arrow', [0.2, 0.3], [0.5, 0.5], ...
% %            'Color', 'r', 'LineWidth', 1.2, 'HeadStyle', 'vback2', ...
% %            'LineStyle', '--');
% 
% saveas(gcf, 'Cross_Slice_y_0p45_compare', 'epsc');
% 
% 
% % y=0.5
% figure;
% plot(Data_r.R_slice1_51, Data_r.U_slice1_51, 'LineWidth', 1.5);
% hold on;
% plot(Data_485.R_slice1_51, Data_485.U_slice1_51, 'LineWidth', 1.5);
% plot(Data_cg.R_slice1_51, Data_cg.U_slice1_51, 'LineWidth', 1.5);
% plot(Slice_Res_ldg2, Slice_Ues_ldg2, 'LineWidth', 1.5);
% hold off;
% legend('Reference', 'RDFM-TPFA', 'RDFM-CG', 'RDFM-LDG', 'Location', 'northeast');
% %局部图
% zoom_x_start = 0.749;  % x1
% zoom_x_end = 0.753;    % x2
% zoom_y_start = 0.494;  % y1
% zoom_y_end = 0.501;    % y2
% 
% % hold on;
% % rectangle('Position', [zoom_x_start, zoom_y_start, zoom_x_end - zoom_x_start, zoom_y_end - zoom_y_start],'EdgeColor', 'r', 'LineWidth', 3.5, 'LineStyle', '-');
% % hold off;
% inset_ax = axes('Position', [0.15, 0.15, 0.3, 0.3]);
% box on;
% grid off;
% set(gca, 'XTick',[],'YTick',[]);
% hold on;
% set(inset_ax, 'FontSize', 8);
% plot(inset_ax, Data_r.R_slice1_51, Data_r.U_slice1_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_485.R_slice1_51, Data_485.U_slice1_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_cg.R_slice1_51, Data_cg.U_slice1_51, 'LineWidth', 1.5);
% plot(inset_ax, Slice_Res_ldg2, Slice_Ues_ldg2, 'LineWidth', 1.5);
% hold(inset_ax, 'off');
% xlim(inset_ax, [zoom_x_start, zoom_x_end]);
% ylim(inset_ax, [zoom_y_start, zoom_y_end]);
% % annotation('arrow', [0.2, 0.3], [0.5, 0.5], ...
% %            'Color', 'r', 'LineWidth', 1.2, 'HeadStyle', 'vback2', ...
% %            'LineStyle', '--');
% 
% saveas(gcf, 'Cross_Slice_y_0p5_compare', 'epsc');
% 
% % x=0.3
% figure;
% plot(Data_r.R_slice2_51, Data_r.U_slice2_51, 'LineWidth', 1.5);
% hold on;
% plot(Data_485.R_slice2_51, Data_485.U_slice2_51, 'LineWidth', 1.5);
% plot(Data_cg.R_slice2_51, Data_cg.U_slice2_51, 'LineWidth', 1.5);
% plot(Slice_Res_ldg3, Slice_Ues_ldg3, 'LineWidth', 1.5);
% hold off;
% legend('Reference', 'RDFM-TPFA', 'RDFM-CG', 'RDFM-LDG', 'Location', 'northeast');
% %局部图
% zoom_x_start = 0.492;  % x1
% zoom_x_end = 0.508;    % x2
% zoom_y_start = 0.5;  % y1
% zoom_y_end = 0.503;    % y2
% 
% % hold on;
% % rectangle('Position', [zoom_x_start, zoom_y_start, zoom_x_end - zoom_x_start, zoom_y_end - zoom_y_start],'EdgeColor', 'r', 'LineWidth', 3.5, 'LineStyle', '-');
% % hold off;
% inset_ax = axes('Position', [0.15, 0.15, 0.3, 0.3]);
% box on;
% grid off;
% set(gca, 'XTick',[],'YTick',[]);
% hold on;
% set(inset_ax, 'FontSize', 8);
% plot(inset_ax, Data_r.R_slice2_51, Data_r.U_slice2_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_485.R_slice2_51, Data_485.U_slice2_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_cg.R_slice2_51, Data_cg.U_slice2_51, 'LineWidth', 1.5);
% plot(inset_ax, Slice_Res_ldg3, Slice_Ues_ldg3, 'LineWidth', 1.5);
% hold(inset_ax, 'off');
% xlim(inset_ax, [zoom_x_start, zoom_x_end]);
% ylim(inset_ax, [zoom_y_start, zoom_y_end]);
% % annotation('arrow', [0.2, 0.3], [0.5, 0.5], ...
% %            'Color', 'r', 'LineWidth', 1.2, 'HeadStyle', 'vback2', ...
% %            'LineStyle', '--');
% 
% saveas(gcf, 'Cross_Slice_x_0p3_compare', 'epsc');
% 
% % x=0.4
% figure;
% plot(Data_r.R_slice3_51, Data_r.U_slice3_51, 'LineWidth', 1.5);
% hold on;
% plot(Data_485.R_slice3_51, Data_485.U_slice3_51, 'LineWidth', 1.5);
% plot(Data_cg.R_slice3_51, Data_cg.U_slice3_51, 'LineWidth', 1.5);
% plot(Slice_Res_ldg4, Slice_Ues_ldg4, 'LineWidth', 1.5);
% hold off;
% legend('Reference', 'RDFM-TPFA', 'RDFM-CG', 'RDFM-LDG', 'Location', 'northeast');
% %局部图
% zoom_x_start = 0.492;  % x1
% zoom_x_end = 0.508;    % x2
% zoom_y_start = 0.5;  % y1
% zoom_y_end = 0.501;    % y2
% 
% % hold on;
% % rectangle('Position', [zoom_x_start, zoom_y_start, zoom_x_end - zoom_x_start, zoom_y_end - zoom_y_start],'EdgeColor', 'r', 'LineWidth', 3.5, 'LineStyle', '-');
% % hold off;
% inset_ax = axes('Position', [0.15, 0.15, 0.3, 0.3]);
% box on;
% grid off;
% set(gca, 'XTick',[],'YTick',[]);
% hold on;
% set(inset_ax, 'FontSize', 8);
% plot(inset_ax, Data_r.R_slice3_51, Data_r.U_slice3_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_485.R_slice3_51, Data_485.U_slice3_51, 'LineWidth', 1.5);
% plot(inset_ax, Data_cg.R_slice3_51, Data_cg.U_slice3_51, 'LineWidth', 1.5);
% plot(inset_ax, Slice_Res_ldg4, Slice_Ues_ldg4, 'LineWidth', 1.5);
% hold(inset_ax, 'off');
% xlim(inset_ax, [zoom_x_start, zoom_x_end]);
% ylim(inset_ax, [zoom_y_start, zoom_y_end]);
% % annotation('arrow', [0.2, 0.3], [0.5, 0.5], ...
% %            'Color', 'r', 'LineWidth', 1.2, 'HeadStyle', 'vback2', ...
% %            'LineStyle', '--');
% 
% saveas(gcf, 'Cross_Slice_x_0p4_compare', 'epsc');
% 
% 
% 
