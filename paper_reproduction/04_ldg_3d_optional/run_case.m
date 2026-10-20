clear;close all;
root=fileparts(mfilename('fullpath'));addpath(root);
N=100;
S=load(fullfile(fileparts(root),'02_single_fracture_3d','fracture_parameters.mat'));
result=solve_ldg(N,S.ParametersofFracture(1,3:14));
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'ldg_result.mat'),'result','-v7.3');
t=linspace(0,1,1001)';xyz=[t 1-t 1-t];
p=evaluate_ldg_pressure(result,xyz);
writetable(table(sqrt(3)*t,p,'VariableNames',{'ArcLength','Pressure'}),fullfile(out,'ldg_profile.csv'));
figure;plot(sqrt(3)*t,p,'k-');xlabel('Arc length');ylabel('Dimensionless pressure');box on;
