function plot_realcase
root=fileparts(mfilename('fullpath'));
file=fullfile(root,'results','mass_balance_data.mat');
if ~exist(file,'file'),file=fullfile(root,'saved_results','mass_balance_data.mat');end
D=load(file);G=D.G;R=abs(D.residual);nx=D.nx;ny=D.ny;
out=fullfile(root,'figures');if ~exist(out,'dir'),mkdir(out);end
x=linspace(0,700,nx+1);y=linspace(0,600,ny+1);[X,Y]=ndgrid(x,y);
Z=reshape(R,nx,ny);Z=Z([1:nx nx],[1:ny ny]);
masks={true(G.cells.num,1),~D.fractured,D.fractured};
names={'all_cell_residuals','matrix_cell_residuals','fracture_cell_residuals'};
for k=1:3
    f=figure('Color','w','Position',[80 80 760 610]);
    mask=reshape(masks{k},nx,ny);A=double(mask([1:nx nx],[1:ny ny]));
    surf(X,Y,Z,'EdgeColor','none','FaceColor','flat','AlphaData',A,'FaceAlpha','flat');
    view(2);axis([0 700 0 600]);pbaspect([700 600 1]);box on;grid off;
    set(gca,'XTick',[],'YTick',[],'FontName','Times New Roman','FontSize',13);
    colormap jet;cb=colorbar;ylabel(cb,'|R_i|');caxis([0 max(R(masks{k}))]);
    set(f,'PaperPositionMode','auto');print(f,fullfile(out,[names{k} '.jpg']),'-djpeg','-r300');
    savefig(f,fullfile(out,[names{k} '.fig']));close(f);
end
f=figure('Color','w','Position',[80 80 760 610]);
plot(1:G.cells.num,R,'.','MarkerSize',4);box on;grid off;
xlabel('Cell index');ylabel('|R_i|');xlim([1 G.cells.num]);
set(gca,'FontName','Times New Roman','FontSize',13);
set(f,'PaperPositionMode','auto');print(f,fullfile(out,'all_cell_residuals_by_index.jpg'),'-djpeg','-r300');
savefig(f,fullfile(out,'all_cell_residuals_by_index.fig'));close(f);
end
