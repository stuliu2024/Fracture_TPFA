function [Tf,report] = axis_aligned_fracture_halves(G,vertices,kappa)
N=G.cartDims;h=1./N;
nc=G.cells.num; nf=size(vertices,1);
cn=rldecode((1:nc)',diff(G.cells.facePos));cf=G.cells.faces(:,1);
fc=G.faces.centroids(cf,:);cc=G.cells.centroids;
delta=fc-cc(cn,:);[~,faceaxis]=max(abs(delta),[],2);
Tf=zeros(numel(cf),1);report.area=zeros(nf,1);report.on_grid_face=false(nf,1);
tol=1e-12;
for f=1:nf
    v=reshape(vertices(f,:),3,[])';lo=min(v,[],1);hi=max(v,[],1);
    normal=find(hi-lo<tol);assert(numel(normal)==1,'Only axis-aligned rectangles supported.');
    tang=setdiff(1:3,normal);pos=lo(normal);
    ix=pos/h(normal);
    if abs(ix-round(ix))<tol && pos>tol && pos<1-tol
        weight=0.5*double(abs(cc(:,normal)-pos)<h(normal)/2+tol);
        report.on_grid_face(f)=true;
    else
        weight=double(pos>cc(:,normal)-h(normal)/2-tol & pos<cc(:,normal)+h(normal)/2+tol);
    end
    lens=zeros(nc,2);
    for j=1:2
        a=tang(j);
        lens(:,j)=max(0,min(cc(:,a)+h(a)/2,hi(a))-max(cc(:,a)-h(a)/2,lo(a)));
    end
    active=weight>0 & all(lens>tol,2);
    report.area(f)=sum(weight(active).*prod(lens(active,:),2));
    expected=prod(hi(tang)-lo(tang));
    assert(abs(report.area(f)-expected)<1e-10,'Fracture area counted incorrectly.');
    for j=1:2
        a=tang(j);other=3-j;
        inside=fc(:,a)>lo(a)+tol & fc(:,a)<hi(a)-tol;
        if lo(a)<tol,inside=inside | abs(fc(:,a))<tol;end
        if hi(a)>1-tol,inside=inside | abs(fc(:,a)-1)<tol;end
        take=faceaxis==a & active(cn) & inside;
        ids=find(take);cells=cn(ids);
        Tf(ids)=Tf(ids)+kappa*weight(cells).*lens(cells,other)/(h(a)/2);
    end
end
assert(all(isfinite(Tf)) && all(Tf>=0));
report.total_area=sum(report.area);
report.nonzero_half_faces=nnz(Tf);
end
