function solve_cg(n)
root=fileparts(mfilename('fullpath'));
addpath(fileparts(root));setup_paths;
xa=0; xb=1;
yc=0; yd=1;
ze=0; zf=1;

nx=n; ny=n; nz=n;

hx=(xb-xa)/nx; hy=(yd-yc)/ny; hz=(zf-ze)/nz;

x=linspace(xa,xb,nx+1);
y=linspace(yc,yd,ny+1);
z=linspace(ze,zf,nz+1);
[xx,yy,zz]=meshgrid(x,y,z); %横向编号
xx=permute(xx,[2,1,3]);yy=permute(yy,[2,1,3]);zz=permute(zz,[2,1,3]);
node=[xx(:),yy(:),zz(:)];

elem = zeros(8,nx*ny*nz);
cellidx = 1:nx*ny*nz;
[i, j, k] = ind2sub([nx ny nz],cellidx);  % index of cells in subscript form
elem(1,cellidx) = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) +i;
elem(2,cellidx) = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) +(i+1);
elem(3,cellidx) = (k-1)*(nx+1)*(ny+1) + ((j+1)-1)*(nx+1) +(i+1);
elem(4,cellidx) = (k-1)*(nx+1)*(ny+1) + ((j+1)-1)*(nx+1) +i;
elem(5,cellidx) = ((k+1)-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) +i;
elem(6,cellidx) = ((k+1)-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) +(i+1);
elem(7,cellidx) = ((k+1)-1)*(nx+1)*(ny+1) + ((j+1)-1)*(nx+1) +(i+1);
elem(8,cellidx) = ((k+1)-1)*(nx+1)*(ny+1) + ((j+1)-1)*(nx+1) +i;

Cell_N = nx; Cell_M = ny; Cell_L = nz;
xmin = xa; xmax = xb; ymin = yc; ymax = yd; zmin = ze; zmax = zf;

[coordinates,elements,EtoEmap] = CubicMesh( xmin,xmax,ymin,ymax,zmin,zmax,Cell_N,Cell_M,Cell_L );% nonperiodic/periodic boundary condition

Tvolumn =  hx*hy*hz; Farea = [hx*hy;hx*hy;hy*hz;hy*hz;hz*hx;hz*hx];
Jacobimat = [1/hx,0,0 ; 0,1/hy,0 ; 0,0,1/hz];%Jacobimat=[Dr/Dx,Dr/Dy,Dr/Dz ; Ds/Dx,Ds/Dy,Ds/Dz ; Dt/Dx,Dt/Dy,Dt/Dz]
FacetoLocalNodes = [1,2,3,4; 5,6,7,8; 1,4,5,8; 2,3,6,7; 1,2,5,6; 3,4,7,8];
neumann = zeros( 2*( nx*ny + ny*nz + nz*nx) , 4) ; neumann_count = 0;
boundary_faces = zeros( 2*( nx*ny + ny*nz + nz*nx) , 4) ; boundary_count = 0;

NumberofFractures = 1;
ParametersofFracture = zeros(NumberofFractures, 18);
ParametersofFracture(:,1) = 1e-3; % width of fractures
ParametersofFracture(:,2) = 1e8; % permeability of fractures
shared=load(fullfile(root,'fracture_parameters.mat'));
ParametersofFracture(:,3:14)=shared.ParametersofFracture(:,3:14); % identical TPFA geometry

HighOrderDegree = 1 ; % bilinear Lagrange basis
HighOrderNp = (HighOrderDegree+1)^3; % degree of freedom in each element for high order polynomial
Gauss_num = 3 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1);
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
[quad_lamda1,quad_lamda2,quad_lamda3]=meshgrid(Gauss_x,Gauss_x,Gauss_x);
quad_lamda1 = quad_lamda1(:);
quad_lamda2 = quad_lamda2(:);
quad_lamda3 = quad_lamda3(:);
quad_w = kron(Gauss_w',Gauss_w*Gauss_w');
quad_w = quad_w(:);
basisP = zeros( HighOrderNp , size(quad_w,1) ) ; % basis function's value in quadrature point
GradrbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
GradsbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
GradtbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
boundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradrboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradsboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradtboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
basisP(1,:) = (1-quad_lamda1).*(1-quad_lamda2).*(1-quad_lamda3) ;
basisP(2,:) =    quad_lamda1 .*(1-quad_lamda2).*(1-quad_lamda3) ;
basisP(3,:) =    quad_lamda1 .*   quad_lamda2 .*(1-quad_lamda3) ;
basisP(4,:) = (1-quad_lamda1).*   quad_lamda2 .*(1-quad_lamda3) ;
basisP(5,:) = (1-quad_lamda1).*(1-quad_lamda2).*quad_lamda3 ;
basisP(6,:) =    quad_lamda1 .*(1-quad_lamda2).*quad_lamda3 ;
basisP(7,:) =    quad_lamda1 .*   quad_lamda2 .*quad_lamda3 ;
basisP(8,:) = (1-quad_lamda1).*   quad_lamda2 .*quad_lamda3 ;

GradrbasisP(1,:) = -(1-quad_lamda2).*(1-quad_lamda3);
GradrbasisP(2,:) =  (1-quad_lamda2).*(1-quad_lamda3);
GradrbasisP(3,:) =     quad_lamda2 .*(1-quad_lamda3) ;
GradrbasisP(4,:) =    -quad_lamda2 .*(1-quad_lamda3) ;
GradrbasisP(5,:) = -(1-quad_lamda2).*quad_lamda3;
GradrbasisP(6,:) =  (1-quad_lamda2).*quad_lamda3;
GradrbasisP(7,:) =     quad_lamda2 .*quad_lamda3 ;
GradrbasisP(8,:) =    -quad_lamda2 .*quad_lamda3;

GradsbasisP(1,:) = -(1-quad_lamda1).*(1-quad_lamda3);
GradsbasisP(2,:) =    -quad_lamda1.*(1-quad_lamda3);
GradsbasisP(3,:) =     quad_lamda1.*(1-quad_lamda3);
GradsbasisP(4,:) =   (1-quad_lamda1).*(1-quad_lamda3);
GradsbasisP(5,:) = -(1-quad_lamda1).*quad_lamda3;
GradsbasisP(6,:) =    -quad_lamda1.*quad_lamda3;
GradsbasisP(7,:) =     quad_lamda1.*quad_lamda3;
GradsbasisP(8,:) =  (1-quad_lamda1).*quad_lamda3;

GradtbasisP(1,:) = -(1-quad_lamda1).*(1-quad_lamda2);
GradtbasisP(2,:) = -quad_lamda1.*(1-quad_lamda2);
GradtbasisP(3,:) = -quad_lamda1.*quad_lamda2;
GradtbasisP(4,:) = -(1-quad_lamda1).*quad_lamda2;
GradtbasisP(5,:) = (1-quad_lamda1).*(1-quad_lamda2);
GradtbasisP(6,:) = quad_lamda1.*(1-quad_lamda2);
GradtbasisP(7,:) = quad_lamda1.*quad_lamda2;
GradtbasisP(8,:) = (1-quad_lamda1).*quad_lamda2;

face_w = Gauss_w*Gauss_w'; face_w=face_w(:);
quad_face_r = zeros(6,Gauss_num^2);
quad_face_s = zeros(6,Gauss_num^2);
quad_face_t = zeros(6,Gauss_num^2);

for ii = 1 : 6
    if ii==1
        [face_r,face_s,face_t] = meshgrid(Gauss_x,Gauss_x,0);
    elseif ii==2
        [face_r,face_s,face_t] = meshgrid(Gauss_x,Gauss_x,1);
    elseif ii==3
        [face_r,face_s,face_t] = meshgrid(0,Gauss_x,Gauss_x);
    elseif ii==4
        [face_r,face_s,face_t] = meshgrid(1,Gauss_x,Gauss_x);
    elseif ii==5
        [face_r,face_s,face_t] = meshgrid(Gauss_x,0,Gauss_x);
    elseif ii==6
        [face_r,face_s,face_t] = meshgrid(Gauss_x,1,Gauss_x);
    end
    face_r=face_r(:);face_s=face_s(:);face_t=face_t(:);
    quad_face_r(ii,:) = face_r;
    quad_face_s(ii,:) = face_s;
    quad_face_t(ii,:) = face_t;
    boundaryP(1,:,ii) = (1-face_r).*(1-face_s).*(1-face_t);
    boundaryP(2,:,ii) = face_r.*(1-face_s).*(1-face_t);
    boundaryP(3,:,ii) = face_r.*face_s.*(1-face_t);
    boundaryP(4,:,ii) = (1-face_r).*face_s.*(1-face_t);
    boundaryP(5,:,ii) = (1-face_r).*(1-face_s).*face_t;
    boundaryP(6,:,ii) = face_r.*(1-face_s).*face_t;
    boundaryP(7,:,ii) = face_r.*face_s.*face_t;
    boundaryP(8,:,ii) = (1-face_r).*face_s.*face_t;

    GradrboundaryP(1,:,ii) = -(1-face_s).*(1-face_t);
    GradrboundaryP(2,:,ii) = (1-face_s).*(1-face_t);
    GradrboundaryP(3,:,ii) = face_s.*(1-face_t);
    GradrboundaryP(4,:,ii) = -face_s.*(1-face_t);
    GradrboundaryP(5,:,ii) = -(1-face_s).*face_t;
    GradrboundaryP(6,:,ii) = (1-face_s).*face_t;
    GradrboundaryP(7,:,ii) = face_s.*face_t;
    GradrboundaryP(8,:,ii) = -face_s.*face_t;

    GradsboundaryP(1,:,ii) = -(1-face_r).*(1-face_t);
    GradsboundaryP(2,:,ii) = -face_r.*(1-face_t);
    GradsboundaryP(3,:,ii) = face_r.*(1-face_t);
    GradsboundaryP(4,:,ii) = (1-face_r).*(1-face_t);
    GradsboundaryP(5,:,ii) = -(1-face_r).*face_t;
    GradsboundaryP(6,:,ii) = -face_r.*face_t;
    GradsboundaryP(7,:,ii) = face_r.*face_t;
    GradsboundaryP(8,:,ii) = (1-face_r).*face_t;

    GradtboundaryP(1,:,ii) = -(1-face_r).*(1-face_s);
    GradtboundaryP(2,:,ii) = -face_r.*(1-face_s);
    GradtboundaryP(3,:,ii) = -face_r.*face_s;
    GradtboundaryP(4,:,ii) = -(1-face_r).*face_s;
    GradtboundaryP(5,:,ii) = (1-face_r).*(1-face_s);
    GradtboundaryP(6,:,ii) = face_r.*(1-face_s);
    GradtboundaryP(7,:,ii) = face_r.*face_s;
    GradtboundaryP(8,:,ii) = (1-face_r).*face_s;
end

[ParametersofFracture,FracturesPath,FracturesArea,basisPonFracture, ...
    GradxbasisPonFracture,GradybasisPonFracture,GradzbasisPonFracture,GradsigmabasisPonFracture] = ...
    SetFractureBasisData3d(coordinates,elements,ParametersofFracture,Gauss_x);

b = zeros(size(coordinates,1),1);% right hand side of the scheme
u = zeros(size(coordinates,1),1);% solution vector
IK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);
JK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);
VK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);

for jj = 1 : size( elements,2 ) % assembly K & b
    xyz = [coordinates(elements(1,jj),1)+hx*quad_lamda1,coordinates(elements(1,jj),2)+hy*quad_lamda2,coordinates(elements(1,jj),3)+hz*quad_lamda3];
    stima = Tvolumn * ( ...
        (GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1)+GradtbasisP*Jacobimat(3,1))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1)+GradtbasisP*Jacobimat(3,1))' ...
        + ...
        (GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2)+GradtbasisP*Jacobimat(3,2))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2)+GradtbasisP*Jacobimat(3,2))'  ...
        + ...
        (GradrbasisP*Jacobimat(1,3)+GradsbasisP*Jacobimat(2,3)+GradtbasisP*Jacobimat(3,3))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,3)+GradsbasisP*Jacobimat(2,3)+GradtbasisP*Jacobimat(3,3))' ) ;
    [J_e,I_e] = meshgrid(elements(:,jj)); I_index = I_e(:); J_index = J_e(:);
    IK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=I_index;
    JK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=J_index;
    VK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=stima(:);
    b(elements(:,jj)) = b(elements(:,jj)) + Tvolumn*basisP*(f(xyz).*quad_w) ;
    for ii = 1 : 6  % the ii-th face of jj-th element
        if ii==1
            outer_vec = [0;0;-1];
        elseif ii==2
            outer_vec = [0;0;1];
        elseif ii==3
            outer_vec = [-1;0;0];
        elseif ii==4
            outer_vec = [1;0;0];
        elseif ii==5
            outer_vec = [0;-1;0];
        else
            outer_vec = [0;1;0];
        end
        E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj
        if E == -1
            boundary_count = boundary_count + 1; boundary_faces(boundary_count,:) = elements(FacetoLocalNodes(ii,:),jj) ;
            if abs(outer_vec(1))<0.5 %  Neumann boundary condition
                neumann_count = neumann_count + 1; neumann(neumann_count,:) = elements(FacetoLocalNodes(ii,:),jj) ;
                Boundary_xyz = ...
                    [coordinates(elements(1,jj),1)+hx*quad_face_r(ii,:)',coordinates(elements(1,jj),2)+hy*quad_face_s(ii,:)',coordinates(elements(1,jj),3)+hz*quad_face_t(ii,:)'];
                b(elements(FacetoLocalNodes(ii,:),jj)) = b(elements(FacetoLocalNodes(ii,:),jj)) + ...
                    Farea(ii)*boundaryP(FacetoLocalNodes(ii,:),:,ii)*(g(Boundary_xyz,outer_vec).*face_w);
            end
        end
    end
end

for k = 1 : NumberofFractures
    for ell = 1 : ParametersofFracture(k,18)
        CurrentIndex = sum(ParametersofFracture(1:k-1,18)) + ell;
        jj = FracturesPath( CurrentIndex );
        [J_e,I_e] = meshgrid(elements(:,jj)); I_index = I_e(:); J_index = J_e(:);
        stima = ParametersofFracture(k,1)*ParametersofFracture(k,2)*FracturesArea( CurrentIndex ) * ( ...
            GradxbasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradxbasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
            + GradybasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradybasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
            + GradzbasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradzbasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
            - GradsigmabasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradsigmabasisPonFracture((CurrentIndex-1)*8+[1:8],:)');
        IK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = I_index;
        JK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = J_index;
        VK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = stima(:);
    end
end

K = sparse(IK,JK,VK,size(coordinates,1),size(coordinates,1));
dirichlet = setdiff(boundary_faces, neumann, 'rows') ;
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
u(unique(dirichlet)) = p(coordinates(unique(dirichlet),:));
b = b - K * u;

clear IK JK VK;
A=K(FreeNodes,FreeNodes);rhs=b(FreeNodes);clear K;
tic;
scale=1./sqrt(diag(A));S=spdiags(scale,0,numel(FreeNodes),numel(FreeNodes));
B=S*A*S;bb=scale.*rhs;
L=ichol((B+B')/2,struct('type','ict','droptol',1e-3,'diagcomp',1e-3));
x0=(1-mod(FreeNodes(:)-1,n+1)/n)./scale;
[y,flag,relres,iterations]=pcg(B,bb,1e-9,3000,L,L',x0);
u(FreeNodes)=scale.*y;solve_seconds=toc;
relative_residual=norm(A*u(FreeNodes)-rhs)/norm(rhs);
if flag~=0,error('CG linear solve did not converge (flag %d).',flag);end
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'cg_result.mat'),'u','n','ParametersofFracture','relative_residual','solve_seconds','flag','iterations','-v7.3');
fprintf('CG: %d free DOFs; relative residual %.6e\n',numel(FreeNodes),relative_residual);
end
function val = km(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
val = ones(size(x));
end
function val = f(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
val = zeros(size(x));
end
function val = p(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
val = 1-x;
end
function val = g(xyz,n)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
val = zeros(size(x));
end

function [coordinates,elements,EtoEmap] = CubicMesh( xmin,xmax,ymin,ymax,zmin,zmax,N,M,L )
XX=linspace(xmin,xmax,N+1); YY=linspace(ymin,ymax,M+1); ZZ=linspace(zmin,zmax,L+1);
[X,Y,Z]=meshgrid(XX,YY,ZZ);      % 纵线横线交错生成方格网节点的x,y,z坐标，分别存入X,Y,Z矩阵
X = permute(X,[2,1,3]); Y = permute(Y,[2,1,3]); Z = permute(Z,[2,1,3]);
coordinates=zeros((N+1)*(M+1)*(L+1),3); elements=zeros(8,N*M*L); EtoEmap=zeros(6,N*M*L);
coordinates(:,1)=X(:);coordinates(:,2)=Y(:);coordinates(:,3)=Z(:);
for K=1:L
    for I=1:M
        for J=1:N
            num = (K-1)*N*M + (I-1)*N + J;
            start=(K-1)*(N+1)*(M+1) + (I-1)*(N+1) + J;
            elements(:,num)=[start;start+1;start+N+2;start+N+1;...
                start+(N+1)*(M+1);start+1+(N+1)*(M+1);start+N+2+(N+1)*(M+1);start+N+1+(N+1)*(M+1)];
            EtoEmap(:,num)=[num-N*M;num+N*M;num-1;num+1;num-N;num+N];
        end
    end
end
EtoEmap(1,1:N*M)=-1;
EtoEmap(2,end-N*M+1:end)=-1;
EtoEmap(3,1:N:end)=-1;
EtoEmap(4,N:N:end)=-1;
EtoEmap(5,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1))=-1 ;
EtoEmap(6,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1))=-1;
isperiodic = 0;
if isperiodic == 1
EtoEmap(1,1:N*M)=(1:N*M)+(L-1)*N*M;
EtoEmap(2,end-N*M+1:end)=1:N*M;
EtoEmap(3,1:N:end)=(1:N:N*M*L)+N-1;
EtoEmap(4,N:N:end)=(1:N:N*M*L);
EtoEmap(5,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1))=...
    reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1);
EtoEmap(6,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1))=...
    reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1);
end

end

function [ParametersofFracture,FracturesPath,FracturesArea,basisPonFracture, ...
    GradxbasisPonFracture,GradybasisPonFracture,GradzbasisPonFracture,GradsigmabasisPonFracture] =...
    SetFractureBasisData3d(coordinates,elements,ParametersofFracture,quad_x)
quad_num = size(quad_x,1);
xmax = max(coordinates(:,1)); xmin = min(coordinates(:,1));
ymax = max(coordinates(:,2)); ymin = min(coordinates(:,2));
zmax = max(coordinates(:,3)); zmin = min(coordinates(:,3));
Cell_N = (elements(4,1)-elements(1,1)-1) ;
Cell_M = (elements(5,1)-elements(1,1)) / (Cell_N+1) - 1 ;
Cell_L = size(elements,2) / (Cell_N*Cell_M) ;
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ; hz = ( zmax - zmin ) / Cell_L ;
NumberofFractures = size(ParametersofFracture,1);
for k=1:NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2)); %叉积，类似于外法向量
    sigma = sigma/norm(sigma,2); %归一化
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3))<0.1  %垂直于XOY平面
        [Ix,Jy,~,~,~]=SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2),quad_x);
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        ParametersofFracture(k,18) = size(Kz,1)*size(Ix,1);
    elseif abs(sigma(1))<0.1  %垂直于YOZ平面
        [Jy,Kz,sy,tz,FracturesLength] = SetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
            coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3),quad_x);
        Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
        ParametersofFracture(k,18) = size(Ix,1)*size(Jy,1);
    elseif abs(sigma(2))<0.1  %垂直于ZOX平面
        [Kz,Ix,tz,rx,FracturesLength] = SetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1),quad_x);
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        ParametersofFracture(k,18) = size(Jy,1)*size(Ix,1);
    else
        error('wrong')
    end
end

FracturesPath = zeros(sum(ParametersofFracture(:,18)),1);
FracturesArea = zeros(sum(ParametersofFracture(:,18)),1);
basisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradxbasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradybasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradzbasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradsigmabasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);

for k=1:NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2));
    sigma = sigma/norm(sigma,2);
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3)<0.1) %垂直于XOY
        [Ix,Jy,rx,sy,FracturesLength] = SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2),quad_x);
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        I = repmat(Ix,size(Kz,1),1);
        J = repmat(Jy,size(Kz,1),1);
        K = kron(Kz,ones(size(Ix,1),1));
        r = repmat(rx,size(Kz,1),quad_num);
        s = repmat(sy,size(Kz,1),quad_num);
        if size(Kz,1)==1
            tz = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (coor_vertex3(3)-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            t = repmat(tz,size(Ix,1),1);
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA = FracturesA *(coor_vertex3(3)-coor_vertex2(3));
        elseif size(Kz,1)==2
            tz1 = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (Kz(1)*hz-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            tz2 = 0 + (coor_vertex3(3)-(Kz(end)-1)*hz)/hz * kron(quad_x',ones(1,quad_num));
            t = [repmat(tz1,size(Ix,1),1);repmat(tz2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
        else
            tz1 = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (Kz(1)*hz-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            tz2 = 0 + (coor_vertex3(3)-(Kz(end)-1)*hz)/hz * kron(quad_x',ones(1,quad_num));
            t = [repmat(tz1,size(Ix,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Ix,1)*(size(Kz,1)-2),1);repmat(tz2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
            FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) = FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) *hz;
        end
        elseif abs(sigma(1))<0.1 % perpendicular to YOZ
            [Jy,Kz,sy,tz,FracturesLength] = SetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
                coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3),quad_x);
            Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
            J = repmat(Jy,size(Ix,1),1);
            K = repmat(Kz,size(Ix,1),1);
            I = kron(Ix,ones(size(Jy,1),1));
            s = repmat(sy,size(Ix,1),quad_num);  %积分点
            t = repmat(tz,size(Ix,1),quad_num);  %积分点
                    if size(Ix,1)==1
            rx = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (coor_vertex3(1)-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            r = repmat(rx,size(Jy,1),1);
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA = FracturesA *(coor_vertex3(1)-coor_vertex2(1));
            elseif size(Ix,1)==2
            rx1 = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (Ix(1)*hx-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            rx2 = 0 + (coor_vertex3(1)-(Ix(end)-1)*hx)/hx * kron(quad_x',ones(1,quad_num));
            r = [repmat(rx1,size(Jy,1),1);repmat(rx2,size(Jy,1),1)];
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
                    else
            rx1 = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (Ix(1)*hx-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            rx2 = 0 + (coor_vertex3(1)-(Ix(end)-1)*hx)/hx * kron(quad_x',ones(1,quad_num));
            r = [repmat(rx1,size(Jy,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Jy,1)*(size(Ix,1)-2),1);repmat(rx2,size(Jy,1),1)];
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
            FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) = FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) *hx;
                    end
         elseif abs(sigma(2))<0.1 % perpendicular to ZOX
        [Kz,Ix,tz,rx,FracturesLength] = SetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1),quad_x);
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        K = repmat(Kz,size(Jy,1),1);
        I = repmat(Ix,size(Jy,1),1);
        J = kron(Jy,ones(size(Kz,1),1));
        t = repmat(tz,size(Jy,1),quad_num);
        r = repmat(rx,size(Jy,1),quad_num);
        if size(Jy,1)==1
            sy = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (coor_vertex3(2)-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            s = repmat(sy,size(Ix,1),1);
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA = FracturesA *(coor_vertex3(2)-coor_vertex2(2));
        elseif size(Jy,1)==2
            sy1 = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (Jy(1)*hy-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            sy2 = 0 + (coor_vertex3(2)-(Jy(end)-1)*hy)/hy * kron(quad_x',ones(1,quad_num));
            s = [repmat(sy1,size(Ix,1),1);repmat(sy2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Jy(1)*hy-coor_vertex2(2));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(2)-(Jy(end)-1)*hy);
        else
            sy1 = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (Jy(1)*hy-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            sy2 = 0 + (coor_vertex3(2)-(Jy(end)-1)*hy)/hy * kron(quad_x',ones(1,quad_num));
            s = [repmat(sy1,size(Ix,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Ix,1)*(size(Jy,1)-2),1);repmat(sy2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Jy(1)*hy-coor_vertex2(2));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(2)-(Jy(end)-1)*hy);
            FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) = FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) *hy;
        end
    else
        error('wrong')
    end

    FracturesPath(sum(ParametersofFracture(1:k-1,18))+1:sum(ParametersofFracture(1:k,18)),1) = ...
    (K-1)*Cell_N*Cell_M + (J-1)*Cell_N + I;
FracturesArea(sum(ParametersofFracture(1:k-1,18))+1:sum(ParametersofFracture(1:k,18)),1) = ...
    FracturesA;

basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s).*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s).*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s.*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s.*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s).*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s).*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s.*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s.*t;

GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-s).*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
     (1-s).*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    s.*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    -s.*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-s).*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    (1-s).*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    s.*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    -s.*t/hx;

GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*t/hy;

GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s/hz;

GradsigmabasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:) = ...
    sigma(1)*GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:)...
   +sigma(2)*GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:)...
   +sigma(3)*GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:) ;

end

function [Ix,Jy,r,s,FracturesLength] = SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,xa,ya,xb,yb,quad_x)
[coordinates2,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );
StartEndElements = zeros(1,2);
StartEndElements(1) = ( ceil((ya-ymin)/hy)-1 )*Cell_N + ceil((xa-xmin)/hx);
StartEndElements(2) = ( ceil((yb-ymin)/hy)-1 )*Cell_N + ceil((xb-xmin)/hx);
CurrentElement = StartEndElements(1);
NumElePassed = 1;
pp = 0;
while ( CurrentElement ~= StartEndElements(2) )
[~, ~, ii_int] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
   ,[xa,xb],[ya,yb]); edge_int = ii_int(:,1);
qq = setdiff(edge_int,pp);
NextElement = EtoEmap(qq,CurrentElement);
pp = find(EtoEmap(:,NextElement) == CurrentElement) ;
CurrentElement = NextElement;
NumElePassed = NumElePassed + 1;
end
Ix = zeros(NumElePassed,1);
Jy = zeros(NumElePassed,1);
r = zeros(NumElePassed,size(quad_x,1));
s = zeros(NumElePassed,size(quad_x,1));
FracturesLength = zeros(NumElePassed,1);
CurrentElement = StartEndElements(1); CurrentIndex = 1;
Ix(CurrentIndex) = mod(CurrentElement-1,Cell_N)+1;
Jy(CurrentIndex) = ceil(CurrentElement/Cell_N);
pp = 0; pplamda = ([xa;ya]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
while ( CurrentElement ~= StartEndElements(2) )
    [x_int, y_int, ii_int] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
       ,[xa,xb],[ya,yb]); edge_int = ii_int(:,1);
    qq = setdiff(edge_int,pp); qqlamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
      if length(edge_int) == 1
            FracturesLength(CurrentIndex) = ...
                sqrt( (x_int-xa)^2 + (y_int-ya)^2 );
      elseif length(edge_int) == 2
            FracturesLength(CurrentIndex) = sqrt( diff(x_int)^2 + diff(y_int)^2 );
      else
          error('intersection points is not one or two');
      end
    r(CurrentIndex,:) = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x] ;
    s(CurrentIndex,:) = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x] ;
    NextElement = EtoEmap(qq,CurrentElement);
    pp = find(EtoEmap(:,NextElement) == CurrentElement); pplamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates2(elements4(1,NextElement),:)')./[hx;hy];
    CurrentElement = NextElement; CurrentIndex = CurrentIndex + 1;
Ix(CurrentIndex) = mod(CurrentElement-1,Cell_N)+1;
Jy(CurrentIndex) = ceil(CurrentElement/Cell_N);
end
[x_int, y_int, ~] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
   ,[xa,xb],[ya,yb]);
qqlamda = ([xb;yb]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
FracturesLength(CurrentIndex) = ...
       sqrt( (x_int-xb)^2 + (y_int-yb)^2 );
r(CurrentIndex,:) = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x] ;
s(CurrentIndex,:) = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x] ;

end

end

function [coordinates,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,N,M )
XX=linspace(xmin,xmax,N+1); YY=linspace(ymin,ymax,M+1);
[X,Y]=meshgrid(XX,YY);      % 纵线横线交错生成方格网节点的x、y坐标，分别存入X、Y矩阵
coordinates=zeros((N+1)*(M+1),2); elements4=zeros(4,N*M); EtoEmap=zeros(4,N*M);
X=X'; Y=Y';
coordinates(:,1)=X(:);coordinates(:,2)=Y(:);
for I=1:M
    for J=1:N
        num = (I-1)*N+J;
        start=(I-1)*(N+1)+J;
        elements4(:,num)=[start;start+1;start+N+2;start+N+1];
        EtoEmap(:,num)=[num-N;num+1;num+N;num-1];
    end
end
EtoEmap(1,1:N)=-1;
EtoEmap(3,(M-1)*N+1:M*N)=-1;
EtoEmap(4,1:N:end)=-1;
EtoEmap(2,N:N:end)=-1;
EtoEmap(1,1:N)=(M-1)*N+1:M*N;
EtoEmap(3,(M-1)*N+1:M*N)=1:N;
EtoEmap(4,1:N:end)=N:N:M*N;
EtoEmap(2,N:N:end)=1:N:M*N;
end
