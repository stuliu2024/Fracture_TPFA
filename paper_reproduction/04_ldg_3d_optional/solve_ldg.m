function result=solve_ldg(N,vertices)
root=fileparts(mfilename('fullpath'));addpath(fileparts(root));setup_paths;

total_clock=tic; total_cpu=cputime;
h=1/N; nd=2*N; n=nd^3; mass=(h/2)^3;
alpha=5/h^3; % Same penalty as dimensionless Cross_test/my_rdfm_12 and Regular/my_rdfm_14.
[Dx,Px,d,e]=ldg_1d(N,true);
[Dn,Pn]=ldg_1d(N,false);
M=speye(nd)*(h/2);
C=kron(speye(N),sparse([1 1;-1/sqrt(3) 1/sqrt(3)]/sqrt(2)));
MC=M*C; Mn=C'*MC;
Bx=kron(MC,kron(MC,Dx*C)); By=kron(MC,kron(Dn*C,MC)); Bz=kron(Dn*C,kron(MC,MC));
Px=C'*Px*C; Pn=C'*Pn*C;
P=alpha*(kron(Mn,kron(Mn,Px))+kron(Mn,kron(Pn,Mn))+kron(Pn,kron(Mn,Mn)));
one=zeros(nd,1); one(1:2:end)=h/sqrt(2);
D=kron(one,kron(one,d)); Ed=alpha*kron(C'*one,kron(C'*one,C'*e));

A=(Bx'*Bx+By'*By+Bz'*Bz)/mass+P;
rhs=-Bx'*D/mass+Ed;
clear P Dx Dn Px Pn;

    data=single_fracture_quadrature(N,vertices);
    nf=numel(data.cells); nq=numel(data.weights);
    rows=zeros(nf*nq*8,1); cols=rows; vals=rows;
    [aa,bb,cc]=ndgrid(0:1,0:1,0:1);
    for k=1:nf
        [ix,iy,iz]=ind2sub([N N N],data.cells(k));
        ids=(2*ix-1+aa(:))+nd*(2*iy-2+bb(:))+nd^2*(2*iz-2+cc(:));
        [rr,jj]=ndgrid((k-1)*nq+(1:nq),ids);
        v=bsxfun(@times,data.phi((k-1)*8+(1:8),:)',sqrt(data.areas(k)*data.weights));
        q=(k-1)*nq*8+(1:nq*8); rows(q)=rr(:);cols(q)=jj(:);vals(q)=v(:);
    end
    Q=sparse(rows,cols,vals,nf*nq,n); clear rows cols vals;
    normal=data.normal(:)/norm(data.normal);
    [~,minaxis]=min(abs(normal)); v=zeros(3,1);v(minaxis)=1;
    t1=cross(normal,v);t1=t1/norm(t1);t2=cross(normal,t1);
    kappa=1e-3*1e8;
    for t=[t1 t2]
        H=(t(1)*(Q*Bx)+t(2)*(Q*By)+t(3)*(Q*Bz))/mass;
        rhs=rhs-kappa*(H'*(t(1)*(Q*D)/mass));
        A=A+kappa*(H'*H);
    end
    fracture_area=sum(data.areas);
    clear Q H data;
clear Bx By Bz D Ed;
[pressure,flag,relres,iterations]=agmg(A,rhs,10,1e-10,3000,0);
relative_residual=norm(A*pressure-rhs)/norm(rhs);
if flag~=0 || relative_residual>1e-7,error('LDG linear solve failed: flag=%d, residual=%g.',flag,relative_residual);end
result=struct('pressure',pressure,'N',N,'vertices',vertices,'pressure_dofs',n, ...
    'seconds',toc(total_clock),'relative_residual',relative_residual,'iterations',iterations, ...
    'fracture_area',fracture_area,'alpha',alpha);

end

function [B,P,d,e]=ldg_1d(N,isdir)
l=[1;-sqrt(3)]/sqrt(2);r=[1;sqrt(3)]/sqrt(2);
S=[0 0;sqrt(3) 0]; B=kron(speye(N),sparse(S));P=sparse(2*N,2*N);
d=zeros(2*N,1);e=d;
for k=1:N-1
    I=2*k-1:2*k;J=2*k+1:2*k+2;
    B(I,I)=B(I,I)-r*r'/2; B(I,J)=B(I,J)-r*l'/2;
    B(J,J)=B(J,J)+l*l'/2; B(J,I)=B(J,I)+l*r'/2;
    P(I,I)=P(I,I)+r*r';P(J,J)=P(J,J)+l*l';
    P(I,J)=P(I,J)-r*l';P(J,I)=P(J,I)-l*r';
end
L=1:2;R=2*N-1:2*N;
if isdir
    P(L,L)=P(L,L)+l*l';P(R,R)=P(R,R)+r*r';
    d(L)=l;e(L)=l; % p_D=1 on x=0, p_D=0 on x=1.
else
    B(L,L)=B(L,L)+l*l';B(R,R)=B(R,R)-r*r';
end
end
