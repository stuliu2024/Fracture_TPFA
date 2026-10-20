function plot_geometry
root=fileparts(mfilename('fullpath'));
c=regular_case21;
out=fullfile(root,'figures');if ~exist(out,'dir'),mkdir(out);end
fig=figure('Color','w','Visible','off','Position',[100 100 700 610]);hold on;
colors=[.30 .55 .75;.76 .53 .30;.48 .65 .48];
for f=1:9
    v=reshape(c.vertices(f,:),3,[])';normal=find(max(v)-min(v)<1e-12);
    opacity=.22;if f>3,opacity=.38;end;if f>6,opacity=.6;end
    ph=patch(v(:,1),v(:,2),v(:,3),colors(normal,:),'FaceAlpha',opacity,'EdgeColor',[.15 .15 .15],'LineWidth',.9);
    if f==1,fracture_handle=ph;end
end
for a=1:3
    other=setdiff(1:3,a);
    for b=[0 1]
        for d=[0 1]
            v=zeros(2,3);v(:,a)=[0;1];v(:,other(1))=b;v(:,other(2))=d;
            plot3(v(:,1),v(:,2),v(:,3),'Color',[.25 .25 .25],'LineWidth',.8);
        end
    end
end
for a=1:3
    other=setdiff(1:3,a);q=[0 0;1 0;1 1;0 1];
    v=zeros(4,3);v(:,other)=q*.25;
    hin=patch(v(:,1),v(:,2),v(:,3),[.1 .4 .8],'FaceAlpha',.75,'EdgeColor','none');
    v=ones(4,3);v(:,other)=.875+q*.125;
    hout=patch(v(:,1),v(:,2),v(:,3),[.65 .18 .55],'FaceAlpha',.85,'EdgeColor','none');
end
hl=plot3([0 1],[0 1],[0 1],'k--','LineWidth',1.3);
axis equal;axis([0 1 0 1 0 1]);view(38,23);box on;grid off;
axis off;
text(-.05,0,.15,'$q_N=-1$','Interpreter','latex','FontSize',18, ...
    'HorizontalAlignment','right','Clipping','off');
text(1.04,1,1.02,sprintf('$h=%g$',c.outlet_head),'Interpreter','latex','FontSize',18, ...
    'HorizontalAlignment','left','Clipping','off');
add_direction_axes(fig,gca);
save_jpg(fig,out,'fracture_geometry');close(fig);
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
