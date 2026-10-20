function plot_si_profiles
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'figures');if ~exist(out,'dir'),mkdir(out);end
methods={'UiB__TPFA','UiB__MPFA','UiB__MVEM','UiB__RT0', ...
    'USTUTT__MPFA','USTUTT__TPFA_Circ','LANL__MFD', ...
    'UNICE_UNIGE__VAG_Cont','UNICE_UNIGE__HFV_Cont', ...
    'UNICE_UNIGE__VAG_Disc','UNICE_UNIGE__HFV_Disc', ...
    'ETHZ_USI__FEM_LM','UNICAMP__Hybrid_Hdiv', ...
    'UNIL_USI__FE_AMR_AFC','DTU__FEM_COMSOL'};
labels={'UiB-TPFA','UiB-MPFA','UiB-MVEM','UiB-RT0', ...
    'USTUTT-MPFA','USTUTT-TPFA Circ','LANL-MFD', ...
    'UNICE UNIGE-VAG Cont','UNICE UNIGE-HFV Cont', ...
    'UNICE UNIGE-VAG Disc','UNICE UNIGE-HFV Disc', ...
    'ETHZ USI-FEM LM','UNICAMP-Hybrid Hdiv', ...
    'UNIL USI-FE AMR AFC','DTU-FEM COMSOL'};
palette=[31 119 180;255 127 14;44 160 44;214 39 40;148 103 189; ...
    140 86 75;227 119 194;127 127 127;188 189 34;23 190 207]/255;
ci=[2 3 1 4 5 6 3 9 7 10 8 4 6 5 8];
ls={'-','-','-',':','-','-','--','-','-',':',':','--','--','--','--'};

ref=read_profile(fullfile(root,'references','reference_USTUTT_MPFA_fine.csv'));
Ns=[9 17 35];
counts={'500','4,000','32,000'};
actualCounts={'729','4,913','42,875'};
fig=figure('Color','w','Visible','off','Units','pixels','Position',[50 50 1200 620]);
axs=gobjects(3,1);
for level=0:2
    ax=axes('Parent',fig,'Position',[.07+.30*level .14 .30 .59]);
    axs(level+1)=ax;hold(ax,'on');h=gobjects(17,1);
    for j=1:numel(methods)
        d=read_profile(fullfile(root,'references',sprintf('level%d',level),[methods{j} '.csv']));
        h(j)=plot(ax,d(:,1),d(:,2),'Color',palette(ci(j),:), ...
            'LineStyle',ls{j},'LineWidth',.85);
    end
    h(16)=plot(ax,ref(:,1),ref(:,2),'k-','LineWidth',1.1);
    N=Ns(level+1);

    r=load_regular_result(N);c=((1:N)-.5)/N;t=linspace(0,1,1001)';
    F=griddedInterpolant({c,c,c},reshape(r.head,N,N,N),'linear','nearest');
    head=F(t,t,t);head(end)=1;
    h(17)=plot(ax,sqrt(3)*t,head,'b--','LineWidth',1.1,'Marker','none');
    xlabel(ax,'Arc length (m)');
    title(ax,{['Published: approximately ' counts{level+1} ' cells'], ...
        sprintf('RDFM-TPFA: %d^3 (%s cells)',N,actualCounts{level+1})},'FontWeight','normal','FontSize',11);
    xlim(ax,[0 sqrt(3)]);ylim(ax,[.5 2.75]);
    set(ax,'FontName','Times New Roman','FontSize',12,'Box','on', ...
        'XTick',0:.5:1.5,'YTick',.5:.25:2.75,'GridColor',[.65 .65 .65], ...
        'GridAlpha',.45,'LineWidth',.7);
    grid(ax,'on');
    if level==0
        ylabel(ax,'Hydraulic head (m)');
        legendHandles=h;
    else
        set(ax,'YTickLabel',[]);
    end
    assert(numel(findobj(ax,'Type','line'))==17);
end
allLabels=[labels {'Reference (USTUTT-MPFA, fine)','RDFM-TPFA'}];
leg=legend(axs(1),legendHandles,allLabels,'NumColumns',4,'Interpreter','none', ...
    'FontName','Times New Roman','FontSize',10,'Box','on');
set(leg,'Units','normalized','Position',[.06 .825 .92 .16]);
for j=1:3,set(axs(j),'Position',[.07+.30*(j-1) .14 .30 .59]);end
set(fig,'PaperUnits','inches','PaperPositionMode','manual', ...
    'PaperPosition',[0 0 12 6.2],'PaperSize',[12 6.2]);
name='head_profiles_9_17_35_three_panels';
print(fig,fullfile(out,[name '.jpg']),'-djpeg','-r300');
savefig(fig,fullfile(out,[name '.fig']));
info=imfinfo(fullfile(out,[name '.jpg']));
assert(info.Width==3600 && info.Height==1860);
close(fig);
end
