function make_efficiency_table
root=fileparts(mfilename('fullpath'));
Nx=[70;140;280;560;1120];Ny=[60;120;240;480;960];
TPFA_DOF=Nx.*Ny;CG_DOF=(Nx-1).*(Ny+1);LDG_DOF=4*Nx.*Ny;
% Recorded times from the original plot_t.m, not a new timing run.
TPFA_s=[.8152407;1.4756511;2.9178926;7.3777868;24.876059];
CG_s=[.8677412;1.6477227;3.9877352;11.6995743;43.5365765];
LDG_s=[1.7260515;5.7950906;23.0910633;148.2262497;422.8074398];
T=table(Nx,Ny,TPFA_DOF,TPFA_s,CG_DOF,CG_s,LDG_DOF,LDG_s);
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
writetable(T,fullfile(out,'cpu_dof_records.csv'));
fid=fopen(fullfile(out,'cpu_dof_rows.tex'),'w');
for k=1:height(T)
    fprintf(fid,'$%d \\times %d$ & %d & %.1f & %d & %.1f & %d & %.1f \\\\\n', ...
        Nx(k),Ny(k),TPFA_DOF(k),TPFA_s(k),CG_DOF(k),CG_s(k),LDG_DOF(k),LDG_s(k));
end
fclose(fid);disp(T);
end
