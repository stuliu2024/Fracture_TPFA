
%SingleFracture3D
clc,clear
close all
mrstModule add incomp linearsolvers
gravity reset off;

%% Initialize a grid
N=40; 
Nx = [N N N];
G = cartGrid(Nx,[1 1 1],'cellnodes',true);
G = computeGeometry(G);
% plotGrid(G);

%% 
westFaces = find(G.faces.centroids(:,1) == 0);
bc = addBC([], westFaces, 'pressure',        ...
           repmat(1, [numel(westFaces), 1]));      
xMax = max(G.faces.centroids(:,1));
eastFaces = find(G.faces.centroids(:,1) == xMax);
bc = addBC(bc, eastFaces, 'pressure',        ...
           zeros(numel(eastFaces), 1));

state0 = initResSol(G,0,1);
perm = 1;
perm = repmat(perm,G.cells.num,1);
rock = struct('perm',perm);
fluid = initSingleFluid('mu' , 1, ...
                        'rho', 1);
T    = computeTrans(G, rock);
[T_fracs,posfracs] = my_computeTrans_fracs_3D(G);
T = T + sparse(posfracs,ones(size(posfracs)),T_fracs,size(T,1),size(T,2));
tic
spparms('spumoni',2) 
state = incompTPFA(state0, G, T,fluid, 'bc',bc);
toc
p              = state.pressure;

% 3d slices
xmin = 0; xmax = 1; ymin = 0; ymax = 1; zmin = 0; zmax =1;
Cell_N = G.cartDims(1); Cell_M = G.cartDims(2); Cell_L = G.cartDims(3);
[XX,YY,ZZ] = meshgrid( linspace(xmin,xmax,Cell_N), linspace(ymin,ymax,Cell_M), linspace(zmin,zmax,Cell_L) );
VV = reshape(p,Cell_N,Cell_M,Cell_L);
VV = permute(VV,[2,1,3]);
xslice = [];   
yslice = [];
zslice = [0.2,0.4,0.6,0.8];
figure;slice(XX,YY,ZZ,VV,xslice,yslice,zslice);
colormap('jet');colorbar;
shading interp
xlim([xmin,xmax]); ylim([ymin,ymax]); zlim([zmin,zmax]);

%y_slice
xslice = [];   
yslice = [0.2,0.4,0.6,0.8];
zslice = [];
figure;slice(XX,YY,ZZ,VV,xslice,yslice,zslice);
colormap('jet');colorbar;
shading interp
xlim([xmin,xmax]); ylim([ymin,ymax]); zlim([zmin,zmax]);

%x_slice
xslice = [0.2,0.4,0.6,0.8];   
yslice = [];
zslice = [];
figure;slice(XX,YY,ZZ,VV,xslice,yslice,zslice);
colormap('jet');colorbar;
shading interp
xlim([xmin,xmax]); ylim([ymin,ymax]); zlim([zmin,zmax]);


%center
xslice = [0.5];   
yslice = [0.5];
zslice = [0.5];
figure;slice(XX,YY,ZZ,VV,xslice,yslice,zslice);
colormap('jet');colorbar;view([-15,25]);
shading interp
xlim([xmin,xmax]); ylim([ymin,ymax]); zlim([zmin,zmax]);

%% function
function [T_fracs,posfracs] = my_computeTrans_fracs_3D(G)
xmin = 0; xmax = 1; ymin = 0; ymax = 1; zmin = 0; zmax =1;
Cell_N = G.cartDims(1); Cell_M = G.cartDims(2); Cell_L = G.cartDims(3);
coordinates = G.nodes.coords;
elements = G.cellNodes(:,[1 2 4 3 5 6 8 7])';

%%
% ParametersofFracture: width,permeability,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sigma1,sigma2,sigma3,...
%                       # of elements passed by this fracture. 

NumberofFractures = 1; 
ParametersofFracture = zeros(NumberofFractures, 18); 
ParametersofFracture(:,1) = 1e-3; % width of fractures
ParametersofFracture(:,2) = 1e8; % permeability of fractures
ParametersofFracture(:,3:14) = [0.1,0.85,0.8,  0.1,0.15,0.2,  0.9,0.15,0.2,  0.9,0.85,0.8]...
                              + (rand(NumberofFractures,12)-0)*1e-4; % position of fractures
[ParametersofFracture,FracturesPath,FracturesArea] = ...
    mySetFractureBasisData3d(coordinates,elements,ParametersofFracture);
% plot the fractures and grids
figure; hold on;
fill3(ParametersofFracture(:,3:3:14)',ParametersofFracture(:,4:3:14)',ParametersofFracture(:,5:3:14)',[0.2392, 0.1490, 0.6588]);
xlim([xmin,xmax]); ylim([ymin,ymax]);zlim([zmin,zmax]); view([-45,20]);grid on; box on;
xlabel('x', 'FontSize', 12)
ylabel('y', 'FontSize', 12)
zlabel('z', 'FontSize', 12)
hold off;
posfracs = (FracturesPath(:,1)'-1)*6+[1:6]';
posfracs = posfracs(:);
facefracs = G.cells.faces(posfracs ,1);
cellNo = rldecode(1:G.cells.num,diff(G.cells.facePos),2).';
cf = G.cells.faces(:,1);
sgn = 2*(cellNo == G.faces.neighbors(cf,1)) - 1; 
N = bsxfun(@times,sgn,G.faces.normals(cf,:));
N_fracs = N(posfracs,:)./sqrt(N(posfracs,1).^2+N(posfracs,2).^2+N(posfracs,3).^2);
sigma_fracs = ParametersofFracture(1:NumberofFractures,15:17);
sst = sigma_fracs'*sigma_fracs;  %sigma*sigmaT
vvt = eye(size(G.cartDims,2)) - sst; % I - sigma*sigmaT 
eplisonkf = ParametersofFracture(:,1)*ParametersofFracture(:,2);
%1./abs(v*n);
v_vel = sqrt([vvt(1,1) vvt(2,2) vvt(3,3)]);   %abs([v1 v2 v3])
v_vel = repmat(v_vel,length(posfracs),1);
kfvn = 1./abs(v_vel(:,1).*N_fracs(:,1)+v_vel(:,2).*N_fracs(:,2)+v_vel(:,3).*N_fracs(:,3));

kfvvt = kfvn.*eplisonkf.*[vvt(1,1) vvt(1,2) vvt(1,3) vvt(2,2) vvt(2,3) vvt(3,3)];
K_Fracs = kfvvt(:, [1, 2, 3, 2, 4, 5, 3, 5, 6]); 
r =      [1, 1, 1, 2, 2, 2, 3, 3, 3] ;
c =      [1, 2, 3, 1, 2, 3, 1, 2, 3] ;

C = G.cells.centroids;
C = G.faces.centroids(G.cells.faces(:,1),:) - C(cellNo,:);
C_fracs = C(posfracs,:);
sigma_fracs1 = sigma_fracs([2 1 3]);
sigma_fracs1 = repmat(sigma_fracs1,2,1); 
sigma_fracs1 = reshape(sigma_fracs1,[],1);
sigma_fracs1 = repmat(sigma_fracs1,length(FracturesPath),1);
area_fracs1 = rldecode(FracturesArea,6,1);
w_fracs = area_fracs1.* sigma_fracs1 ./G.faces.areas(facefracs,:);
w_fracs = abs(w_fracs);
C_fracs = C_fracs.*w_fracs;
% Compute T = C'*K*N / C'*C.
T_fracs = zeros(size(posfracs,1),1);
for k=1:size(r,2)
    T_fracs = T_fracs + (C_fracs(:,r(k)).*K_Fracs(:,k).*N_fracs(:,c(k)));
end
T_fracs = T_fracs./sum(C_fracs.*C_fracs,2);

end

%% 
function [coordinates,elements,EtoEmap] = CubicMesh( xmin,xmax,ymin,ymax,zmin,zmax,N,M,L )
% time : 2020.2.11 - 2020.2.11
% author : xuziyao
%   本函数目前仅适用于立方体区域的均匀立方体剖分.
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
% periodic boundary condition :
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

function [ParametersofFracture,FracturesPath,FracturesArea] =...
    mySetFractureBasisData3d(coordinates,elements,ParametersofFracture)
% Global geometric information:
xmax = max(coordinates(:,1)); xmin = min(coordinates(:,1));
ymax = max(coordinates(:,2)); ymin = min(coordinates(:,2));
zmax = max(coordinates(:,3)); zmin = min(coordinates(:,3));
% x direction is divided into N parts ,y direction is divided into M parts, z direction is divided into L parts
Cell_N = (elements(4,1)-elements(1,1)-1) ; 
Cell_M = (elements(5,1)-elements(1,1)) / (Cell_N+1) - 1 ;
Cell_L = size(elements,2) / (Cell_N*Cell_M) ; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ; hz = ( zmax - zmin ) / Cell_L ;
% set fractures data:
NumberofFractures = size(ParametersofFracture,1);
% ParametersofFracture: width,permeability,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sigma1,sigma2,sigma3,...
%                       # of elements passed by this fracture. 
for k = 1 : NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2));
    sigma = sigma/norm(sigma,2);
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3))<0.1 % perpendicular to XOY
        [Ix,~,~] = mySetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2));
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        ParametersofFracture(k,18) = size(Kz,1)*size(Ix,1);
    elseif abs(sigma(1))<0.1 % perpendicular to YOZ
        [Jy,~,~] = mySetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
            coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3));
        Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
        ParametersofFracture(k,18) = size(Ix,1)*size(Jy,1);
    elseif abs(sigma(2))<0.1 % perpendicular to ZOX
        [~,Ix,~] = mySetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1));
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        ParametersofFracture(k,18) = size(Jy,1)*size(Ix,1);
    else
        error('wrong')
    end
end

FracturesPath = zeros(sum(ParametersofFracture(:,18)),1);
FracturesArea = zeros(sum(ParametersofFracture(:,18)),1);

for k = 1 : NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2));
    sigma = sigma/norm(sigma,2);
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3))<0.1 % perpendicular to XOY
        [Ix,Jy,FracturesLength] = mySetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2));
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        I = repmat(Ix,size(Kz,1),1);
        J = repmat(Jy,size(Kz,1),1);
        K = kron(Kz,ones(size(Ix,1),1));
        if size(Kz,1)==1
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA = FracturesA *(coor_vertex3(3)-coor_vertex2(3));
        elseif size(Kz,1)==2
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
        else
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
            FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) = FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) *hz;
        end
    elseif abs(sigma(1))<0.1 % perpendicular to YOZ
        [Jy,Kz,FracturesLength] = mySetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
            coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3));
        Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
        J = repmat(Jy,size(Ix,1),1);
        K = repmat(Kz,size(Ix,1),1);
        I = kron(Ix,ones(size(Jy,1),1));
        if size(Ix,1)==1
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA = FracturesA *(coor_vertex3(1)-coor_vertex2(1));
        elseif size(Ix,1)==2
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
        else
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
            FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) = FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) *hx;
        end
    elseif abs(sigma(2))<0.1 % perpendicular to ZOX
        [Kz,Ix,FracturesLength] = mySetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1));
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        K = repmat(Kz,size(Jy,1),1);
        I = repmat(Ix,size(Jy,1),1);
        J = kron(Jy,ones(size(Kz,1),1));
        if size(Jy,1)==1
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA = FracturesA *(coor_vertex3(2)-coor_vertex2(2));
        elseif size(Jy,1)==2
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Jy(1)*hy-coor_vertex2(2));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(2)-(Jy(end)-1)*hy);
        else
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
end


function [Ix,Jy,FracturesLength] = mySetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,xa,ya,xb,yb)
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
end


%% 函数
function A = rldecode(A, n, dim)
%Decompress run length encoding of array `A` along dimension `dim`.
%
% SYNOPSIS:
%   B = rldecode(A, n, dim)
%   B = rldecode(A, n) % dim assumed to be 1
%
% PARAMETERS:
%   A         - encoded array
%   n         - repetition of each layer along dimension `dim`. `n` can be
%               either a scalar or one repetition number for each layer.
%   dim       - dimension of `A` where run length encoding is done.
%               dim > 0.
%
% RETURNS:
%   B         - uncompressed `A`
%
% EXAMPLE:
%   % 1) Numerical example:
%   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
%   [B,n] = rlencode(A,1)
%   C = rldecode(B,n,1)
%
%   % 2) Retrive 'first' column of G.cells.faces (see grid_structure):
%   G = cartGrid([10, 10, 2]);
%   cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
%   disp(['CellFace nr. 10 belongs to cell nr: ', num2str(cellNo(10))]);
%
% SEE ALSO:
%   `rlencode`

%{This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).%}

if nargin < 3
  dim = 1;
end

assert(dim > 0, 'Third argument DIM must be positive');

if numel(n) == 1,
   n = repmat(n, [size(A, dim), 1]);
end

assert (all( n(:)>=0 ), 'All repeat counts should be nonnegative.');
if nargin < 3,
   assert (numel(n) == size(A, dim), ...
   sprintf(['There should be a repeat count for each value along dimension dim.\n',...
    'The default value of dim is 1. Did you forget to specify dim?']));
else
   assert (numel(n) == size(A, dim), ...
   'There should be a repeat count for each value along dimension dim.');
end

% Take dimension we compress along to be first dimension,
% i.e., swap dimensions 1 and dim.
d      = 1:max(dim, ndims(A));
d([1, dim])   = [dim, 1];
B      = permute(A,d);

r      = n(:)~=0;
B      = reshape(B(r, :), sum(r), []);

% Insert repeated layers and permute back
i      = cumsum([1; double(reshape(n(r), [], 1))]);
j      = zeros(i(end)-1,1);
j(i(1:end-1)) = 1;

szA    = [size(A), ones(1, dim-ndims(A))];
A      = permute(reshape(B(cumsum(j),:), [sum(n(:)), szA(d(2:end))]), d);
end


end
