function plot_regular_network
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'figures');if ~exist(out,'dir'),mkdir(out);end
plot_geometry;
files=dir(fullfile(root,'references','published_level2','*__*.csv'));
curves=cell(numel(files),1);lo=-Inf;hi=Inf;
for j=1:numel(files)
    curves{j}=read_profile(fullfile(files(j).folder,files(j).name));
    lo=max(lo,curves{j}(1,1));hi=min(hi,curves{j}(end,1));
end
arc=linspace(lo,hi,1000)';v=zeros(1000,numel(files));
for j=1:numel(files),v(:,j)=interp1(curves{j}(:,1),curves{j}(:,2),arc,'linear');end
v=sort(v,2);bounds=zeros(1000,2);
for j=1:2
    q=[.1 .9];index=1+(size(v,2)-1)*q(j);a=floor(index);b=ceil(index);
    bounds(:,j)=(b-index)*v(:,a)+(index-a)*v(:,b);
    if a==b,bounds(:,j)=v(:,a);end
end
ref=read_profile(fullfile(root,'references','reference_USTUTT_MPFA_fine.csv'));
fig=figure('Color','w','Position',[80 80 760 610]);hold on;
fill([arc;flipud(arc)],[bounds(:,1);flipud(bounds(:,2))],[.76 .76 .76],'EdgeColor','none');
plot(ref(:,1),ref(:,2),'k-','LineWidth',1.8);
colors=[.85 .32 .12;.12 .58 .28;0 .25 .90];styles={':','-.','--'};
t=linspace(0,1,1001)';distance=sqrt(3)*t;meshes=[20 35 50];heads=zeros(numel(t),3);
for j=1:3
    N=meshes(j);r=load_regular_result(N);c=((1:N)-.5)/N;
    F=griddedInterpolant({c,c,c},reshape(r.head,N,N,N),'linear','nearest');
    heads(:,j)=F(t,t,t);heads(end,j)=1;
    plot(distance,heads(:,j),'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',2);
end
xlabel('Arc length (m)');ylabel('Hydraulic head (m)');xlim([0 sqrt(3)]);ylim([.9 2.45]);box on;grid off;
set(gca,'FontName','Times New Roman','FontSize',13);
legend({'Published methods (10-90%, ~32k cells)','Reference (USTUTT-MPFA, fine)', ...
    'RDFM-TPFA (20 x 20 x 20)','RDFM-TPFA (35 x 35 x 35)','RDFM-TPFA (50 x 50 x 50)'}, ...
    'Location','northeast','FontSize',10);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.6 6.1],'PaperSize',[7.6 6.1]);
print(fig,fullfile(out,'head_profile_mesh_comparison.jpg'),'-djpeg','-r300');
savefig(fig,fullfile(out,'head_profile_mesh_comparison.fig'));close(fig);
writetable(array2table([distance heads],'VariableNames',{'ArcLength_m','N20','N35','N50'}),fullfile(out,'profiles_main.csv'));
end
