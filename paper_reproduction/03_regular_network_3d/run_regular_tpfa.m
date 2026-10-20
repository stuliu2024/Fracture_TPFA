function result = run_regular_tpfa(N)
if nargin<1,N=50;end
root=fileparts(mfilename('fullpath'));addpath(root);
addpath(fileparts(root));setup_paths;
out=fullfile(root,'results',sprintf('N%d',N));if ~exist(out,'dir'),mkdir(out);end
c=regular_case21;
tic;G=cartGrid([N N N],[1 1 1]);G=computeGeometry(G);h=1/N;
centers=G.cells.centroids;low_fraction=zeros(G.cells.num,1);
for j=1:size(c.low_boxes,1)
    b=c.low_boxes(j,:);v=ones(G.cells.num,1);
    for a=1:3
        v=v.*max(0,min(centers(:,a)+h/2,b(2*a))-max(centers(:,a)-h/2,b(2*a-1)))/h;
    end
    low_fraction=low_fraction+v;
end
assert(max(low_fraction)<=1+1e-12);
rock.perm=c.matrix_high+(c.matrix_low-c.matrix_high)*low_fraction;
Tmatrix=computeTrans(G,rock);
[Tf,fracture_report]=axis_aligned_fracture_halves(G,c.vertices,c.tangential_conductivity);
T=Tmatrix+Tf;
bf=find(any(G.faces.neighbors==0,2));cent=G.faces.centroids(bf,:);
[~,ax]=max(abs(G.faces.normals(bf,:)),[],2);
inarea=zeros(numel(bf),1);outarea=inarea;
for a=1:3
    tang=setdiff(1:3,a);
    sel=ax==a & abs(cent(:,a))<1e-12;
    inarea(sel)=prod(max(0,min(cent(sel,tang)+h/2,.25)-max(cent(sel,tang)-h/2,0)),2);
    sel=ax==a & abs(cent(:,a)-1)<1e-12;
    outarea(sel)=prod(max(0,min(cent(sel,tang)+h/2,1)-max(cent(sel,tang)-h/2,.875)),2);
end
assert(abs(sum(inarea)-3*.25^2)<1e-12 && abs(sum(outarea)-3*.125^2)<1e-12);
in=bf(inarea>0);outfaces=bf(outarea>0);
bc=addBC([],in,'flux',c.inlet_flux_density*inarea(inarea>0));
bc=addBC(bc,outfaces,'pressure',ones(numel(outfaces),1)*c.outlet_head);
frac=zeros(G.faces.num,1);frac(outfaces)=outarea(outarea>0)./G.faces.areas(outfaces);
cf=G.cells.faces(:,1);cut=frac(cf)>0;
assert(all(Tf(cut)==0),'Outlet patch intersects a fracture: requires explicit edge clipping.');
T(cut)=T(cut).*frac(cf(cut));
setup_seconds=toc;
fluid=initSingleFluid('mu',1,'rho',1);
state=incompTPFA(initResSol(G,0,1),G,T,fluid,'bc',bc,'LinSolve',@(A,b) agmg(A,b,10,1e-10,3000,0),'MatrixOutput',true);
relative_residual=norm(state.A*state.pressure-state.rhs)/norm(state.rhs);
assert(relative_residual<1e-7,'Linear-system residual exceeds 1e-7.');
solve=struct('true_relres',relative_residual,'passed_residual_check',true);
cn=rldecode((1:G.cells.num)',diff(G.cells.facePos));
sgn=2*(G.faces.neighbors(cf,1)==cn)-1;
residual=accumarray(cn,sgn.*state.flux(cf),[G.cells.num 1]);
bsign=2*(G.faces.neighbors(bf,1)>0)-1;
outward=bsign.*state.flux(bf);
result=struct('N',N,'method','RDFM-TPFA axis-aligned branch','head',state.pressure, ...
    'case_settings',c,'fracture_report',fracture_report,'solve',solve, ...
    'setup_seconds',setup_seconds,'total_seconds',toc,'cell_residual',residual, ...
    'max_abs_mass_residual',max(abs(residual)),'net_boundary_outflow',sum(outward), ...
    'inflow',-sum(outward(inarea>0)),'outflow',sum(outward(outarea>0)), ...
    'matrix_low_fraction',low_fraction,'dofs',G.cells.num);
assert(abs(result.inflow-.1875)<1e-10,'Incorrect prescribed inlet total.');
save(fullfile(out,'rdfm_tpfa_solution.mat'),'result','-v7.3');
fprintf('N=%d DOFs=%d head=[%.10g,%.10g] relres=%.4e\n',N,result.dofs,min(result.head),max(result.head),result.solve.true_relres);
fprintf('AREA=%.12g INFLOW=%.12g OUTFLOW=%.12g MAX_CELL_RESIDUAL=%.4e\n',fracture_report.total_area,result.inflow,result.outflow,result.max_abs_mass_residual);

end
