function plot_single_fracture

root=fileparts(mfilename('fullpath'));
data=fullfile(root,'results');
if ~exist(fullfile(data,'tpfa_result.mat'),'file') || ~exist(fullfile(data,'cg_result.mat'),'file'),data=fullfile(root,'saved_results');end
T=load(fullfile(data,'tpfa_result.mat'));
C=load(fullfile(data,'cg_result.mat'));
S=load(fullfile(root,'fracture_parameters.mat'));
assert(T.N==C.n,'TPFA and CG must use the same mesh.');
assert(isequal(C.ParametersofFracture(:,3:14),S.ParametersofFracture(:,3:14)));
assert(all(isfinite(T.p)) && all(isfinite(C.u)));
out=fullfile(root,'figures');if ~exist(out,'dir'),mkdir(out);end
v=reshape(S.ParametersofFracture(1,3:14),3,[])';
fig=figure('Color','w','Visible','off','Position',[50 50 760 610]);hold on;
patch('Vertices',v,'Faces',[1 2 3 4], ...
    'FaceColor',[.30 .55 .75],'FaceAlpha',.65,'EdgeColor',[.2 .2 .2],'LineWidth',1);
for a=1:3
    other=setdiff(1:3,a);
    for b=[0 1]
        for d=[0 1]
            q=zeros(2,3);q(:,a)=[0;1];q(:,other(1))=b;q(:,other(2))=d;
            plot3(q(:,1),q(:,2),q(:,3),'Color','k','LineWidth',.8);
        end
    end
end
q=[0 0;1 0;1 1;0 1];
patch(zeros(4,1),q(:,1),q(:,2),[.1 .4 .8],'FaceAlpha',.08,'EdgeColor','none');
patch(ones(4,1),q(:,1),q(:,2),[.65 .18 .55],'FaceAlpha',.08,'EdgeColor','none');
plot3([0 1],[1 0],[1 0],'k--','LineWidth',1.3);
axis equal;axis([0 1 0 1 0 1]);view(38,23);axis off;
text(-.05,0,.5,'$p_D=1$','Interpreter','latex','FontSize',18,'HorizontalAlignment','right','Clipping','off');
text(1.04,1,.5,'$p_D=0$','Interpreter','latex','FontSize',18,'HorizontalAlignment','left','Clipping','off');
add_direction_axes(fig,gca);
save_jpg(fig,out,'single_fracture_geometry');close(fig);

t=linspace(0,1,1001)';x=t;y=1-t;z=1-t;s=sqrt(3)*t;
N=T.N;c=((1:N)-.5)/N;
F=griddedInterpolant({c,c,c},reshape(T.p,N,N,N),'linear','nearest');tpfa=F(x,y,z);
tpfa(1)=1;tpfa(end)=0;
c=linspace(0,1,N+1);
F=griddedInterpolant({c,c,c},reshape(C.u,N+1,N+1,N+1),'linear','none');cg=F(x,y,z);
fig=figure('Color','w','Visible','off','Position',[50 50 760 610]);hold on;
plot(s,tpfa,'b-','LineWidth',2);plot(s,cg,'r--','LineWidth',2);
xlabel('Arc length');ylabel('Dimensionless pressure');
set(gca,'FontName','Times New Roman','FontSize',13);
xlim([0 sqrt(3)]);box on;grid off;
legend({'RDFM-TPFA','RDFM-CG'},'Location','northeast','FontSize',12);
save_jpg(fig,out,'single_fracture_pressure_profiles');close(fig);
writetable(table(s,x,y,z,tpfa,cg,'VariableNames',{'ArcLength_m','x_m','y_m','z_m','RDFM_TPFA','RDFM_CG'}),fullfile(out,'pressure_profiles.csv'));

end
function add_direction_axes(fig,mainax)
ax=axes('Parent',fig,'Units','normalized','Position',[.055 .055 .16 .16]);
hold(ax,'on');
for k=1:3
    d=zeros(1,3);d(k)=1;
    quiver3(ax,0,0,0,d(1),d(2),d(3),0,'k','LineWidth',1.2,'MaxHeadSize',.35);
    text(ax,1.2*d(1),1.2*d(2),1.2*d(3),char('x'+k-1), ...
        'FontName','Times New Roman','FontSize',14,'HorizontalAlignment','center');
end
axis(ax,'equal');axis(ax,[-.2 1.35 -.2 1.35 -.2 1.35]);
view(ax,get(mainax,'View'));axis(ax,'off');
end
function save_jpg(fig,out,name)
set(fig,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 7.6 6.1],'PaperSize',[7.6 6.1]);
print(fig,fullfile(out,[name '.jpg']),'-djpeg','-r300');savefig(fig,fullfile(out,[name '.fig']));
info=imfinfo(fullfile(out,[name '.jpg']));assert(info.Width==2280 && info.Height==1830);
end
