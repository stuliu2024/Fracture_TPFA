clc,clear
close all

iii = [1:5];

Err_p = zeros(size(iii,2),1);
Err_f = zeros(size(iii,2),1);
h_step = zeros(size(iii,2),1);

%theta
theta = 1; %THETA \in [-pi/2,pi/2]

for iii = iii

n=20*2^(iii-1)+1;
nx=n; ny=n;
xa=-pi;xb=pi;yc=-pi;yd=pi;
hx=(xb-xa)/nx;
hy=(yd-yc)/ny;
tic

%% 
x=linspace(xa,xb,nx+1); x=x(:); sx = numel(x)-1; 
y=linspace(yc,yd,ny+1); y=y(:); sy = numel(y)-1; 
celldim = [sx,sy];
numC = sx*sy;  
numN = (sx+1)*(sy+1); 
numFX = (sx+1)*sy; numFY = sx*(sy+1);
numF = numFX + numFY; 
% Nodes/Coordinates 
[xCoords, yCoords] = ndgrid(x,y);
coords = [xCoords(:), yCoords(:)]; 
% Generate face-edges
N = reshape(1:numN,[sx+1, sy+1]);
NF1 = reshape(N(1:sx+1,1:sy),1,[]); 
NF2 = reshape(N(1:sx+1,2:sy+1),1,[]);
faceNodesX = reshape([NF1;NF2],[],1); 
NF1 = reshape(N(1:sx,1:sy+1),1,[]);
NF2 = reshape(N(2:sx+1,1:sy+1),1,[]);
faceNodesY = reshape([NF2;NF1],[],1);
faceNodes = [faceNodesX;faceNodesY]; 
% Generate cell-faces
foffset = 0;
FX = reshape(foffset+(1:numFX),sx+1,sy);
foffset = foffset + numFX;
FY = reshape(foffset+(1:numFY),sx,sy+1);
F1 = reshape(FX(1:sx,:),1,[]); %West == 4
F2 = reshape(FX(2:sx+1,:),1,[]); %East == 2
F3 = reshape(FY(:,1:sy),1,[]); %South == 1
F4 = reshape(FY(:,2:sy+1),1,[]); %North == 3
cellFaces = [reshape([F3;F2;F4;F1],[],1),kron(ones([numC,1]),[1,2,3,4]')]; 
% Generate neighbors
C = zeros([sx+2,sy+2]);
C(2:sx+1,2:sy+1) = reshape(1:numC,[sx,sy]);
NX1 = reshape(C(1:sx+1,2:sy+1),[],1);
NX2 = reshape(C(2:sx+2,2:sy+1),[],1);
NY1 = reshape(C(2:sx+1,1:sy+1),[],1);
NY2 = reshape(C(2:sx+1,2:sy+2),[],1);
neighbors = [[NX1,NX2];[NY1,NY2]];
%cells neighbors cells
NX3 = reshape(C(1:sx,2:sy+1),[],1);
NX4 = reshape(C(3:sx+2,2:sy+1),[],1); 
NY3 = reshape(C(2:sx+1,1:sy),[],1);
NY4 = reshape(C(2:sx+1,3:sy+2),[],1);
cell2cell = [NY3,NX4,NY4,NX3];
G.cell2cell = cell2cell;
% Generate cell nodes
k=reshape(N(1:nx,1:ny),[],1);
cNodes =[k k+1 k+(sx+2) k+(sx+1)]; 
% Assemble structure
G.cells = struct('num',numC,'facePos',(1:4:(numC+1)*4)','indexMap',(1:numC)');
G.faces = struct('num',numF,'nodePos',(1:2:(numF+1)*2)','neighbors',neighbors,'tag',zeros(numF,1));
G.nodes = struct('num',numN,'coords',coords);
G.cells.faces = cellFaces; G.faces.nodes = faceNodes; G.cellNodes = cNodes;
G.cartDims = celldim; G.griddim = numel(G.cartDims);
edges = reshape(G.faces.nodes,2,[]).'; 
[n1, n2] = deal(G.nodes.coords(edges(:,1),:), G.nodes.coords(edges(:,2),:)); 
edgeLength = n2 - n1; 
faceAreas = sqrt(sum(edgeLength.^2,2));
faceCentroids = (n1+n2)./2;
faceNormals = [edgeLength(:,2), -edgeLength(:,1)];
numfaces = diff(G.cells.facePos);
c = faceCentroids(G.cells.faces(:,1),:); 
w=1; n=numfaces; 
no = rldecode(1:numel(n),n,2).';
accum = sparse(no,1:numel(no),w); c = accum*[c,ones([size(c,1),1])]; 
w = c(:,end); c= bsxfun(@rdivide,c(:,1:end-1),w);
cCenter = c; cellno = no; 
cellEdges = edges(G.cells.faces(:,1),:); 
r = G.faces.neighbors(G.cells.faces(:,1),2) == cellno; 
cellEdges(r,:) = cellEdges(r,[2,1]); 
cc = cCenter(cellno,:);
a = G.nodes.coords(cellEdges(:,1),:) - cc;
b = G.nodes.coords(cellEdges(:,2),:) - cc;
quadArea = @(a,b)abs(a(:,1).*b(:,2)-a(:,2).*b(:,1)); 
subArea = quadArea(a,b)./2; 
subCentroid = (cCenter(cellno,:) + 2*faceCentroids(G.cells.faces(:,1),:))/3; 
c = subCentroid; w=subArea; n = numfaces;
no = rldecode(1:numel(n),n,2).';
accum = sparse(no,1:numel(no),w); c = accum*[c,ones([size(c,1),1])];
w = c(:,end); c= bsxfun(@rdivide,c(:,1:end-1),w);
cellCentroids = c; cellVolumes = w; 
% Update grid
G.faces.areas = faceAreas;
G.faces.normals = faceNormals;
G.faces.centroids = faceCentroids;
G.cells.volumes = cellVolumes;
G.cells.centroids = cellCentroids;
%% 
if (-pi/4<=theta && theta<=pi/4) 
    lengthfracs = (xb-xa)*sqrt(1+tan(theta)^2);
else
    lengthfracs = (yd-yc)*sqrt(1+1/tan(theta)^2); 
end

zxface = G.faces.centroids(:,1); 
maxzx = max(zxface); minzx = min(zxface);
rightfaces = find(G.faces.centroids(:,1) > (maxzx - 1e-6));
leftfaces = find(G.faces.centroids(:,1) < (minzx + 1e-6));
zyface = G.faces.centroids(:,2); 
maxzy = max(zyface); minzy = min(zyface);
topfaces = find(G.faces.centroids(:,2) > (maxzy - 1e-6));
bottomfaces = find(G.faces.centroids(:,2) < (minzy + 1e-6));
nrf = numel(rightfaces);   %Dirichlet
nlf = numel(leftfaces);    
ntf = numel(topfaces);     
nbf = numel(bottomfaces);  
bc = struct('face',[],'type',{{}},'value',[]);
bc.face = [leftfaces; rightfaces; bottomfaces; topfaces]; 
bc.type = repmat({'pressure'},[1,numel(bc.face)]);
coords_bc = G.faces.centroids(bc.face,:);
value = pp(coords_bc,theta);
bc.value =value;

[nc,nf] = deal(G.cells.num, G.faces.num);
state = struct('pressure',zeros([nc,1]), 'flux',zeros([nf,1]));

perm = 1;
perm = repmat(perm,G.cells.num,1);
rock = struct('perm',perm);

%% 
cellno = rldecode(1:G.cells.num,diff(G.cells.facePos),2).'; 
cellNo = cellno;
C = G.cells.centroids; 
C = G.faces.centroids(G.cells.faces(:,1),:) - C(cellNo,:); 
cf = G.cells.faces(:,1); 
sgn = 2*(cellNo == G.faces.neighbors(cf,1)) - 1; 
N = bsxfun(@times,sgn,G.faces.normals(cf,:));
K = rock.perm * [1, 0, 1]; K = K(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; 
% Compute T = C'*K*N / C'*C.
T = zeros(size(cellNo));
for k=1:size(i,2)
    T = T + (C(:,i(k)).*K(cellNo,k).*N(:,j(k)));
end
T = T./sum(C.*C,2);
is_neg = T<0;
if any(is_neg)
    T(is_neg) = -T(is_neg);
end
%% RDFM
%-------------------------Fractures---------------------------------------%
fractures = [0+0.001*hx^3, 0+0.004*hx^3, theta, lengthfracs-0.03*hx^3, 1,2]; %xc,yc,angle,length,thickness,permeability
% ²¹³äxa,xb,ya,yb,
fractures = [fractures, fractures(:,1)-0.5*fractures(:,4).*cos(fractures(:,3)), fractures(:,1)+0.5*fractures(:,4).*cos(fractures(:,3)),...
                        fractures(:,2)-0.5*fractures(:,4).*sin(fractures(:,3)), fractures(:,2)+0.5*fractures(:,4).*sin(fractures(:,3)),];
figure
surf(xCoords,yCoords,zeros(nx+1,ny+1));
% colormap('white');
axis image;axis on;ax=axis;axis(ax*1.02);view([0,90]);grid off;
hold on
plot(fractures(:,[7,8])',fractures(:,[9,10])',...
   'r-','LineWidth',1.5); hold off; drawnow;
startendfracs = zeros(size(fractures,1),2);
Nodes = G.nodes.coords;
Elems = G.cellNodes;
Cell2cell = G.cell2cell;
for j=1:size(Elems,1)
    startfracselem =inpolygon(fractures(:,7),fractures(:,9),Nodes(Elems(j,[1;2;3;4;1]),1),Nodes(Elems(j,[1;2;3;4;1]),2));
    startendfracs(startfracselem,1) = j;
    endfracselem = inpolygon(fractures(:,8),fractures(:,10),Nodes(Elems(j,[1;2;3;4;1]),1),Nodes(Elems(j,[1;2;3;4;1]),2));
    startendfracs(endfracselem,2) = j;
end
fractures(:,11)=0; 
fractures(:,12)=0; 
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    preciousfracelem=currentfracselem;  
    while(currentfracselem ~= startendfracs(i,2))
        [~,~,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1); 
        neighborfracsindex = Cell2cell(currentfracselem,edgenum); 
        nextfracselem = setdiff(neighborfracsindex,preciousfracelem); 
        preciousfracelem = currentfracselem;
        currentfracselem = nextfracselem;
        fractures(i,11)=fractures(i,11)+1;
        fractures(i,12)=fractures(i,12)+size(edgenum,1);
    end
    [~,~,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1);
    fractures(i,11)=fractures(i,11)+1;
    fractures(i,12)=fractures(i,12)+size(edgenum,1);
end

pathfractures = zeros(sum(fractures(:,11)),1);
lengthfractures = zeros(sum(fractures(:,11)),1);
centrfractures = zeros(sum(fractures(:,11)),2); 
currentelemindex = 0;
Facesfracs = []; 
K_fracs = [];
v_vel = [];
len_fracs = [];
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    currentelemindex = currentelemindex + 1;
    pathfractures(currentelemindex) = currentfracselem;
    previousedgenum=0;
    leftfracs = ([fractures(i,7);fractures(i,9)]);
    v=[cos(fractures(i,3)), sin(fractures(i,3))];
    v_vel = [v_vel;[repmat(fractures(i,5)*fractures(i,6),fractures(i,11).*size(G.cellNodes,2),1),repmat([cos(fractures(i,3)), sin(fractures(i,3))],fractures(i,11).*size(G.cellNodes,2),1)]];
    while(currentfracselem ~= startendfracs(i,2)) 
        [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1); 
        if (size(edgenum,1) == 1)
            lengthfractures(currentelemindex) = sqrt((xi-fractures(i,7))^2+(yi-fractures(i,9))^2);
        elseif (size(edgenum,1) == 2)
            lengthfractures(currentelemindex) = sqrt(diff(xi)^2+diff(yi)^2);
        end
        tmp_len = lengthfractures(currentelemindex);
        len_fracs = [len_fracs;[repmat(tmp_len,size(G.cellNodes,2),1),repmat([sin(fractures(i,3));cos(fractures(i,3))],size(G.cellNodes,2)/2,1)]];
        neighborfracsindex=setdiff(edgenum,previousedgenum);
        rightfracs = ([xi(neighborfracsindex == edgenum);yi(neighborfracsindex == edgenum)]);
        centrfractures(currentelemindex,:) = (rightfracs'+leftfracs')/2;
        Facesfracs = [Facesfracs;[repmat(currentfracselem,numel(edgenum),1),edgenum,xi,yi,repmat(centrfractures(currentelemindex,:),numel(edgenum),1)]]; %#ok<*AGROW>
        K_fracs = [K_fracs; [repmat(fractures(i,5)*fractures(i,6),numel(edgenum),1),repmat(v,numel(edgenum),1)]];
        nextfracselem = Cell2cell(currentfracselem,neighborfracsindex); 
        previousedgenum = find(Cell2cell(nextfracselem,:) == currentfracselem); 
        currentfracselem = nextfracselem; 
        leftfracs = ([xi(neighborfracsindex == edgenum);yi(neighborfracsindex == edgenum)]); 
        currentelemindex = currentelemindex + 1; 
        pathfractures(currentelemindex) = currentfracselem;
    end

    [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1);
    rightfracs = ([fractures(i,8);fractures(i,10)]); 
    centrfractures(currentelemindex,:) = (rightfracs'+leftfracs')/2;
    Facesfracs = [Facesfracs;[repmat(currentfracselem,numel(edgenum),1),edgenum,xi,yi,repmat(centrfractures(currentelemindex,:),numel(edgenum),1)]];
    K_fracs = [K_fracs; [repmat(fractures(i,5)*fractures(i,6),numel(edgenum),1),v]];
    lengthfractures(currentelemindex) = sqrt((xi-fractures(i,8))^2+(yi-fractures(i,10))^2);
    tmp_len = lengthfractures(currentelemindex);
    len_fracs = [len_fracs;[repmat(tmp_len,size(G.cellNodes,2),1),repmat([sin(fractures(i,3));cos(fractures(i,3))],size(G.cellNodes,2)/2,1)]];
    
end

posfracs1 = (pathfractures(:,1)'-1)*size(G.cellNodes,2)+[1:4]'; posfracs1 = posfracs1(:); posfracs2 = posfracs1;
facsfracs1 = G.cells.faces(posfracs1(:),1);    facsfracs2 = G.cells.faces(posfracs2,1);
uN_fracs2 = N(posfracs1,:)./sqrt(N(posfracs1,1).^2+N(posfracs1,2).^2);
v_p_n = v_vel(:,2).*uN_fracs2(:,1)+v_vel(:,3).*uN_fracs2(:,2);
active_edges = find(abs(v_p_n) >1e-15);
posfracs1 = posfracs1(active_edges,:);
facsfracs1 = facsfracs1(active_edges,:);
v_vel = v_vel(active_edges,:);
uN_fracs2 = uN_fracs2(active_edges,:);
kfvn = 1./abs(v_vel(:,2).*uN_fracs2(:,1)+v_vel(:,3).*uN_fracs2(:,2));
kfvvt = kfvn.*v_vel(:,1).* ...
    [v_vel(:,2).*v_vel(:,2),v_vel(:,2).*v_vel(:,3),v_vel(:,3).*v_vel(:,3)];
K_Fracs = kfvvt(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2];
C_fracs = C(posfracs1,:);
N_fracs = uN_fracs2;
w_fracs =len_fracs(:,1).*len_fracs(:,2)./G.faces.areas(facsfracs2,:);
w_fracs = abs(w_fracs);
w_fracs = w_fracs(active_edges,:);
C_fracs = C_fracs.*w_fracs;
% Compute T = C'*K*N / C'*C.
T_fracs = zeros(size(facsfracs1,1),1);
for k=1:size(i,2)
    T_fracs = T_fracs + (C_fracs(:,i(k)).*K_Fracs(:,k).*N_fracs(:,j(k)));
end
T_fracs = T_fracs./sum(C_fracs.*C_fracs,2);
%% 
T = T + sparse(posfracs1,ones(size(posfracs1)),T_fracs,size(T,1),size(T,2));
%%
neighborship = G.faces.neighbors;   nif = size(neighborship,1); 
cf = G.cells.faces(:,1);  ncf = size(cf,1);
nc = G.cells.num;  n = nc; 
ft = 1./accumarray(cf,1./T,[nif,1]);
if ~isempty(bc)
    hh=zeros(nif,1);
    ff=zeros(ncf,1); 
    gg=zeros(nc,1);
    dF = false([G.faces.num, 1]);
    dC = [];
    is_press = strcmpi('pressure',bc.type);
    face = bc.face(is_press);
    dC = bc.value(is_press); 
    map = sparse(double(face),1,1:numel(face)); 
    dF(face) = true; 
    i = dF(G.cells.faces(:,1)); 
    ff(i) = -dC(map(G.cells.faces(i,1))); 
    dC = dC(map(dF)); 
    is_flux = strcmpi('flux',bc.type);
    hh(bc.face(is_flux)) = -bc.value(is_flux);
end
i=all(neighborship ~=0, 2);
j = i(cf) | dF(cf); 
rhs = accumarray(cellNo,-ft(cf).*ff,[n,1]) + gg  + accumarray(cellNo,-hh(cf),[n,1]);
d = zeros(G.cells.num,1);
d = d + accumarray(cellNo(dF(cf)),T(dF(cf)),[nc,1]) +...
    accumarray(reshape(neighborship(i,:),[],1), repmat(ft(i),[2,1]), [nc,1]); 
I = [neighborship(i,1); neighborship(i,2); (1:nc)'];
J = [neighborship(i,2); neighborship(i,1); (1:nc)'];
V = [-ft(i); -ft(i); d];
A = sparse(double(I), double(J), V, nc, nc);
p = A\rhs; 
fpress = accumarray(cf, p(cellNo).*T, [nif,1])./accumarray(cf(:,1),T,[nif,1]);
b = any(G.faces.neighbors==0,2);
fpress(b) = fpress(b) - hh(b)./ft(b);
fpress(dF)=dC;
sgn = 2*(G.faces.neighbors(~i,2)==0)-1;
ni = neighborship(i,:); 
flux = -accumarray(find(i), ft(i).*(p(ni(:,2))-p(ni(:,1))),[nif,1]); 
c = sum(G.faces.neighbors(~i,:),2);
flux(~i) = -sgn.*ft(~i).*(fpress(~i) - p(c)); 
state.pressure(1:nc) = p(1:nc);
state.flux = flux;
state.facePressure = fpress;
p_c = pp(G.cells.centroids,theta);
u_c = pp_g(G.faces.centroids,theta);
ku_c = -([1 0;0 1])*u_c'; ku_c = ku_c';
f_c = ku_c(:,1).*G.faces.normals(:,1) + ku_c(:,2).*G.faces.normals(:,2);
%L2
error_p = sqrt(sum(G.cells.volumes.*(p_c-p).^2))/((max(p_c)-min(p_c))*sum(G.cells.volumes));
error_f = sqrt(sum(G.faces.areas.*(f_c - flux).^2))/((max(f_c)-min(f_c))*sum(G.faces.areas));
Err_p(iii) = error_p;
Err_f(iii) = error_f;
h_step(iii) = hx;
figure;
[x_rect, y_rect] = meshgrid(linspace(xa, xb, nx), linspace(yc, yd, ny));
error = abs(reshape(p_c,nx,ny)' - reshape(p,nx,ny)');
surf(x_rect, y_rect,error);shading interp;
colormap('jet');
view([0,90]);colorbar;
end

h_step_1=h_step;
err_p_L2_1=Err_p;
err_f_L2_1=Err_f;
k1=size(h_step_1,1);

figure;
set(gcf,'Units','normal');
set(gcf,'Position',[0.25,0.25,0.55,0.4]);
myshowrateh1(h_step_1(1:k1),err_p_L2_1(1:k1),1,'-*','|| p-p_h||',h_step_1(1:k1),err_f_L2_1(1:k1),1,'-.','|| f-f_h||');

err_p_L2_1
(log(err_p_L2_1(1))-log(err_p_L2_1(2)))/(log(h_step_1(1))-log(h_step_1(2)))
(log(err_p_L2_1(2))-log(err_p_L2_1(3)))/(log(h_step_1(2))-log(h_step_1(3)))
(log(err_p_L2_1(3))-log(err_p_L2_1(4)))/(log(h_step_1(3))-log(h_step_1(4)))
(log(err_p_L2_1(4))-log(err_p_L2_1(5)))/(log(h_step_1(4))-log(h_step_1(5)))

err_f_L2_1
(log(err_f_L2_1(1))-log(err_f_L2_1(2)))/(log(h_step_1(1))-log(h_step_1(2)))
(log(err_f_L2_1(2))-log(err_f_L2_1(3)))/(log(h_step_1(2))-log(h_step_1(3)))
(log(err_f_L2_1(3))-log(err_f_L2_1(4)))/(log(h_step_1(3))-log(h_step_1(4)))
(log(err_f_L2_1(4))-log(err_f_L2_1(5)))/(log(h_step_1(4))-log(h_step_1(5)))

toc
%% º¯Êý
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


function val = pp(xy,theta)  %dirichlet±ß½ç
xy = ([cos(theta),sin(theta);-sin(theta),cos(theta)]*xy')';
xx = xy(:,1);
yy = xy(:,2);
val = sin(xx).*exp(abs(yy));
end


function val = pp_g(xy,theta)  %grad_p
xy = ([cos(theta),sin(theta);-sin(theta),cos(theta)]*xy')';
xx = xy(:,1);
yy = xy(:,2);
val = [cos(theta).*cos(xx).*exp(abs(yy))-sin(theta).*sin(xx).*exp(abs(yy)), sin(theta).*cos(xx).*exp(abs(yy))+cos(theta).*sin(xx).*exp(abs(yy))];
end

