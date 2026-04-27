%Cross case.
clc,clear
close all
nnn=485;
nx=nnn; ny=nnn;
xa=0;xb=1;yc=0;yd=1;
hx=(xb-xa)/nx;
hy=(yd-yc)/ny;
tic
addpath(genpath(pwd))
%% 网格信息
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
subCentroid = (cCenter(cellno,:) + 2*faceCentroids(G.cells.faces(:,1),:))/3; %子单元中心
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
zface = G.faces.centroids(:,1);
maxz = max(zface); minz = min(zface);
rightfaces = find(G.faces.centroids(:,1) > (maxz - 1e-6));
leftfaces = find(G.faces.centroids(:,1) < (minz + 1e-6));
nrf = numel(rightfaces); %Dirichlet
nlf = numel(leftfaces);  %Dirichlet
bc = struct('face',[],'type',{{}},'value',[]);
bc.face = [leftfaces; rightfaces]; 
bc.type = repmat({'pressure'},[1,numel(bc.face)]);
value = zeros(nrf+nlf,1); value(1:nlf) = 1; 
bc.value =value;
[nc,nf] = deal(G.cells.num, G.faces.num);
state = struct('pressure',zeros([nc,1]), 'flux',zeros([nf,1])); 
perm = 1;
perm = repmat(perm,G.cells.num,1);
rock = struct('perm',perm);
%% T
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
fractures = [0.5-1e-9,0.5-1e-9,0,   0.5-4e-8,1e-03,1e8;
             0.5-1e-9,0.5-1e-9,pi/2,0.5-4e-8,1e-03,1e8;
             ]; %xc,yc,angle,length,thickness,permeability  
%xa,xb,ya,yb,
fractures = [fractures, fractures(:,1)-0.5*fractures(:,4).*cos(fractures(:,3)), fractures(:,1)+0.5*fractures(:,4).*cos(fractures(:,3)),...
                        fractures(:,2)-0.5*fractures(:,4).*sin(fractures(:,3)), fractures(:,2)+0.5*fractures(:,4).*sin(fractures(:,3)),];
figure
surf(xCoords,yCoords,zeros(nx+1,ny+1));
colormap('white');
axis image;axis on;ax=axis;axis(ax);view([0,90]);grid off; axis off;
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
    end
    [~,~,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1);
    fractures(i,11)=fractures(i,11)+1;
end
pathfractures = zeros(sum(fractures(:,11)),1);
currentelemindex = 0;
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    currentelemindex = currentelemindex + 1;
    pathfractures(currentelemindex) = currentfracselem;
     previousedgenum=0;
    while(currentfracselem ~= startendfracs(i,2)) 
        [~,~,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1);
        neighborfracsindex=setdiff(edgenum,previousedgenum);
        nextfracselem = Cell2cell(currentfracselem,neighborfracsindex); 
        previousedgenum = find(Cell2cell(nextfracselem,:) == currentfracselem); 
        currentfracselem = nextfracselem; 
        currentelemindex = currentelemindex + 1; 
        pathfractures(currentelemindex) = currentfracselem;
    end
end
frac_nodes = zeros(2*size(pathfractures,1),2);
node_index = 0; 
currentelemindex = 0;
nodefracs = cell(size(fractures,1),1);
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    currentelemindex = currentelemindex + 1;
    pathfractures(currentelemindex) = currentfracselem;
    previousedgenum=0;
    node_index = node_index + 1;
    frac_nodes(node_index,:) = fractures(i,[7,9]); 
    while(currentfracselem ~= startendfracs(i,2))
        [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1); 
        nodefracs{i}=[nodefracs{i};[xi,yi]];
        node_index = node_index + size(edgenum,1); 
        temp_index = unique([node_index - (size(edgenum,1)-1),node_index]); 
        frac_nodes(temp_index ,:)= [xi,yi];
        neighborfracsindex=setdiff(edgenum,previousedgenum);
        nextfracselem = Cell2cell(currentfracselem,neighborfracsindex);
        previousedgenum = find(Cell2cell(nextfracselem,:) == currentfracselem); 
        currentfracselem = nextfracselem; 
        currentelemindex = currentelemindex + 1; 
        pathfractures(currentelemindex) = currentfracselem;
    end
    [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1); 
    nodefracs{i}=[nodefracs{i};[xi,yi]];
    node_index = node_index + size(edgenum,1);
    temp_index = unique([node_index - (size(edgenum,1)-1),node_index]);
    frac_nodes(temp_index ,:)= [xi,yi];
    node_index = node_index + 1;
    frac_nodes(node_index ,:)= fractures(i,[8,10]); 
end
frac_mid_nodes = squeeze(sum(reshape(frac_nodes, 2, size(frac_nodes,1)/2, size(frac_nodes,2)),1))/2;
frac_mid_nodes = frac_mid_nodes + (rand(size(frac_mid_nodes,1),size(frac_mid_nodes,2)))*1e-10;
fracs_local_num = zeros(size(frac_mid_nodes,1),2);
fracs_local_dist = zeros(size(frac_mid_nodes,1),2);
for k=1:size(fractures,1)
    for j=1:fractures(k,11)
        currentelemindex = sum(fractures(1:k-1,11))+j;
        currentelem = pathfractures(currentelemindex,1);
        local_nodes = Nodes(Elems(currentelem,[1;2;3;4]),:);
        x1 = local_nodes(1,1); y1 = local_nodes(1,2);
        x2 = local_nodes(2,1); y2 = local_nodes(2,2);
        x3 = local_nodes(3,1); y3 = local_nodes(3,2);
        x4 = local_nodes(4,1); y4 = local_nodes(4,2);
        local_frac_nodes = frac_mid_nodes(currentelemindex,:);
        x0 = local_frac_nodes(1,1); y0 = local_frac_nodes(1,2);
        distance(1,:) = abs((y1-y2)*x0-(x1-x2)*y0+x1*y2-x2*y1)/sqrt((y1-y2)^2+(x1-x2)^2);
        distance(3,:) = abs((y3-y4)*x0-(x3-x4)*y0+x3*y4-x4*y3)/sqrt((y3-y4)^2+(x3-x4)^2);
        if ( abs(sin(fractures(k,3))) < abs(sin(pi/2)) )
            local1 = find(distance == min(distance([1 3])));
            dist1 = min(distance([1 3]));
        elseif ( abs(sin(fractures(k,3))) >=  abs(sin(pi/2)) )
            local1 = -1;
            dist1 = -1;
        end
        distance(2,:) = abs((y2-y3)*x0-(x2-x3)*y0+x2*y3-x3*y2)/sqrt((y2-y3)^2+(x2-x3)^2);
        distance(4,:) = abs((y4-y1)*x0-(x4-x1)*y0+x4*y1-x1*y4)/sqrt((y4-y1)^2+(x4-x1)^2);
        if ( abs(cos(fractures(k,3))) < abs(cos(0)) )
            local2 = find(distance == min(distance([2 4])));
            dist2 = min(distance([2 4]));
        elseif ( abs(cos(fractures(k,3))) >=  abs(cos(0)) )
            local2 = -1;
            dist2 = -1;
        end
        fracs_local_num(currentelemindex,:) = [local1 local2];  
        fracs_local_dist(currentelemindex,:) = [dist1 dist2];  
    end
end

fracs_local_num1 = reshape(fracs_local_num',[],1);
fracs_local_dist1 = reshape(fracs_local_dist',[],1);
pathfractures1 = rldecode(pathfractures,2,1);
fracs_local_k = rldecode((1:size(fractures,1))',fractures(1:size(fractures,1),11),1);
fracs_local_k1 = rldecode(fracs_local_k,2,1);
actual_fracs_info  = [pathfractures1,fracs_local_num1,fracs_local_dist1,fracs_local_k1];
mask1 = ((actual_fracs_info(:,2) == -1) | (actual_fracs_info(:,3) == -1));
actual_fracs_info1 = actual_fracs_info(~mask1,:);
fracs_global_mask = (actual_fracs_info1(:,1)-1)*4 + actual_fracs_info1(:,2);
fracs_global_num = G.cells.faces(fracs_global_mask,1);
fracs_global_elem = G.faces.neighbors(fracs_global_num,:);
[n,m] = size(fracs_global_elem);
same_mask  = (actual_fracs_info1(:,1) == fracs_global_elem);
[~, sort_idx] = sort(~same_mask, 2);
linear_indices = sub2ind([n, m], repmat((1:n)',1,m), sort_idx);
fracs_global_elem1 = fracs_global_elem(linear_indices);
mask2 = ((fracs_global_elem1(:,2) == 0));                                                                         %new
virtual_fracs_global_elem1_fix = fracs_global_elem1(~mask2,2);                                                    %new
fracs_global_num_fix = fracs_global_num(~mask2);                                                                  %new
face_test = reshape(G.cells.faces((virtual_fracs_global_elem1_fix-1)*4 + [1:4],1),[],4);                          %new
virtual_local_num1 = arrayfun(@(row)find(face_test(row,:) == fracs_global_num_fix(row),1),1:size(face_test,1))';  %new
virtual_fracs_info1 = [virtual_fracs_global_elem1_fix, virtual_local_num1, actual_fracs_info1(~mask2,4)];         %new
actual_local_nums = [mod(actual_fracs_info1(:,2)-2,4)+1 mod(actual_fracs_info1(:,2),4)+1];
virtual_local_nums = [mod(virtual_fracs_info1(:,2)-2,4)+1 mod(virtual_fracs_info1(:,2),4)+1];
actual_local_nums1 = reshape(actual_local_nums',[],1);
actual_fracs_info2 = [rldecode(actual_fracs_info1(:,1),2,1),actual_local_nums1];
virtual_local_nums1 = reshape(virtual_local_nums',[],1);
virtual_fracs_info2 = [rldecode(virtual_fracs_info1(:,1),2,1),virtual_local_nums1];
lengthfractures = zeros(sum(fractures(:,11)),1);
thetafractures = zeros(sum(fractures(:,11)),2);
eplisonkfractures = zeros(sum(fractures(:,11)),1);
currentelemindex = 0;
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    currentelemindex = currentelemindex + 1;
    
    thetafractures(currentelemindex,1) = cos(fractures(i,3));
    thetafractures(currentelemindex,2) = sin(fractures(i,3));
    eplisonkfractures(currentelemindex,1) = fractures(i,5)*fractures(i,6);
    
    previousedgenum=0;
    while(currentfracselem ~= startendfracs(i,2))
        [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1);
        if (size(edgenum,1) == 1)
            lengthfractures(currentelemindex) = sqrt((xi-fractures(i,7))^2+(yi-fractures(i,9))^2);
        elseif (size(edgenum,1) == 2)
            lengthfractures(currentelemindex) = sqrt(diff(xi)^2+diff(yi)^2);
        end
        neighborfracsindex=setdiff(edgenum,previousedgenum);
        nextfracselem = Cell2cell(currentfracselem,neighborfracsindex); 
        previousedgenum = find(Cell2cell(nextfracselem,:) == currentfracselem);
        currentfracselem = nextfracselem; 
        currentelemindex = currentelemindex + 1;
        
        thetafractures(currentelemindex,1) = cos(fractures(i,3));
        thetafractures(currentelemindex,2) = sin(fractures(i,3));
        eplisonkfractures(currentelemindex,1) = fractures(i,5)*fractures(i,6);
    end
    [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
    lengthfractures(currentelemindex) = sqrt((xi-fractures(i,8))^2+(yi-fractures(i,10))^2);
    
end
actual_fracs_flux_info = [rldecode(lengthfractures,2,1), reshape(thetafractures',[],1), rldecode(eplisonkfractures,2,1), rldecode(thetafractures,2,1)];
actual_fracs_flux_info1 = actual_fracs_flux_info(~mask1,:);
actual_fracs_flux_info2 = [actual_fracs_info2, rldecode(actual_fracs_flux_info1,2,1), rldecode(actual_fracs_info1(:,[3 4]),2,1)];
actual_fracs_pos = (actual_fracs_flux_info2(:,1)-1)*4 + actual_fracs_flux_info2(:,2);
actual_fracs_N = N(actual_fracs_pos,:)./sqrt(N(actual_fracs_pos,1).^2+N(actual_fracs_pos,2).^2);
%1/|v*n|
kfvn = 1./abs(actual_fracs_flux_info2(:,6).*actual_fracs_N(:,1)+actual_fracs_flux_info2(:,7).*actual_fracs_N(:,2));
%eplison*k_f*1/|v*n|*v*vt
kfvvt = kfvn.*actual_fracs_flux_info2(:,5).*...
    [actual_fracs_flux_info2(:,6).*actual_fracs_flux_info2(:,6),...
    actual_fracs_flux_info2(:,6).*actual_fracs_flux_info2(:,7),...
    actual_fracs_flux_info2(:,7).*actual_fracs_flux_info2(:,7)];
K_Fracs = kfvvt(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; 
C_fracs = C(actual_fracs_pos,:);
N_fracs = actual_fracs_N;
actual_fracs_face = G.cells.faces(actual_fracs_pos,1);
w_fracs = actual_fracs_flux_info2(:,3).*actual_fracs_flux_info2(:,4)./G.faces.areas(actual_fracs_face,:);
w_fracs = abs(w_fracs);
C_fracs = C_fracs.*w_fracs;
% Compute T = C'*K*N / C'*C.
T_fracs = zeros(size(actual_fracs_pos,1),1);
for k=1:size(i,2)
    T_fracs = T_fracs + (C_fracs(:,i(k)).*K_Fracs(:,k).*N_fracs(:,j(k)));
end
T_fracs = T_fracs./sum(C_fracs.*C_fracs,2);
virtual_fracs_pos = (virtual_fracs_info2(:,1)-1)*4+virtual_fracs_info2(:,2);
C_fracs_len = sqrt(C(actual_fracs_pos,1).^2+C(actual_fracs_pos,2).^2);
w1 = (C_fracs_len + actual_fracs_flux_info2(:,8))./(2*C_fracs_len);
w2 = 1 - w1;
w2 = w2(~rldecode(mask2,2,1),:);                                       %new
w1=ones(size(w1,1),1);
w2=zeros(size(w2,1),1);
w=[w1;w2];
Fracs_pos_all = [actual_fracs_pos;virtual_fracs_pos];
T_fracs_all = w.*[T_fracs;T_fracs(~rldecode(mask2,2,1),:)];            %new
%% T
T = T + sparse(Fracs_pos_all,ones(size(Fracs_pos_all)),T_fracs_all,size(T,1),size(T,2));
%%
neighborship = G.faces.neighbors;   nif = size(neighborship,1); 
cf = G.cells.faces(:,1);  ncf = size(cf,1); 
nc = G.cells.num;  n = nc; 
ft = 1./accumarray(cf,1./T,[nif,1]); 
if ~isempty(bc)
    hh=zeros(nif,1);  %Neumann
    ff=zeros(ncf,1);  %Dirichlet
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
Cpu_time = toc;
%% postprocessing
figure;
[x_rect, y_rect] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
contourf(x_rect,y_rect,reshape(p,nx,ny)',20,'linestyle','no');
hold on
plot(fractures(:,[7,8])',fractures(:,[9,10])',...
   'r-','LineWidth',1.5); hold off;
colorbar;

figure;
surf(x_rect,y_rect,reshape(p,nx,ny)');
view([30 50]);

dof = nx*ny;
uc = p;
elem = Elems;
node = Nodes;
sparsity = nnz(A)/dof^2*1000;
conditionalnumber = condest(A);

matrix_elements = nx*ny
fracture_elements = size(pathfractures,1)
unique_temp = setdiff(1:nx*ny,pathfractures)';
size(unique_temp,1)+size(unique(pathfractures),1) - nx*ny
Points = G.nodes.coords;
T = G.cellNodes;
P = p;
[x_rect, y_rect] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
X_meshgrid = x_rect';
Y_meshgrid = y_rect';
U_meshgrid = reshape(p,nx,ny);
U_fracs = cell(size(nodefracs,1),1);
fracXAndPCell = cell(size(nodefracs,1),1);
for i=1:size(fractures,1)
    nodefracs{i} = sort([fractures(i,[7 9]);nodefracs{i};fractures(i,[8 10])],1);
    U_fracs{i} = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',nodefracs{i}(:,1),nodefracs{i}(:,2),'linear');
    fracXAndPCell{i} = [nodefracs{i},U_fracs{i}];
end

if exist('z', 'var') && ~isempty(z)
    filename = ['Data_3e', num2str(z)];
else
    filename = ['Data_', num2str(nnn)];
end
save(filename,'uc','dof','elem','node','matrix_elements','fracture_elements','Cpu_time','sparsity','conditionalnumber','Points','T','P','fracXAndPCell');


%slice

% y=0.45
[X_slice0_51, Y_slice0_51] = meshgrid(linspace(xa, xb, nx),0.45);
R_slice0_51 = sqrt((X_slice0_51 - X_slice0_51(1)).^2+(Y_slice0_51-Y_slice0_51(1)).^2);
U_slice0_51 = interp2(x_rect,y_rect,reshape(p,nx,ny)',X_slice0_51,Y_slice0_51,'linear');
save(filename,'X_slice0_51','Y_slice0_51','R_slice0_51','U_slice0_51','-append') 
figure;plot(R_slice0_51,U_slice0_51,'LineWidth',1.5);

% y=0.5
[X_slice1_51, Y_slice1_51] = meshgrid(linspace(xa, xb, nx),0.5);
R_slice1_51 = sqrt((X_slice1_51 - X_slice1_51(1)).^2+(Y_slice1_51-Y_slice1_51(1)).^2);
U_slice1_51 = interp2(x_rect,y_rect,reshape(p,nx,ny)',X_slice1_51,Y_slice1_51,'linear');
save(filename,'X_slice1_51','Y_slice1_51','R_slice1_51','U_slice1_51','-append') 
figure;plot(R_slice1_51,U_slice1_51,'LineWidth',1.5);

% x=0.3
[X_slice2_51, Y_slice2_51] = meshgrid(0.3,linspace(yc, yd, ny));
R_slice2_51 = sqrt((X_slice2_51 - X_slice2_51(1)).^2+(Y_slice2_51-Y_slice2_51(1)).^2);
U_slice2_51 = interp2(x_rect,y_rect,reshape(p,nx,ny)',X_slice2_51,Y_slice2_51,'linear');
save(filename,'X_slice2_51','Y_slice2_51','R_slice2_51','U_slice2_51','-append') 
figure;plot(R_slice2_51,U_slice2_51,'LineWidth',1.5);

% x=0.4
[X_slice3_51, Y_slice3_51] = meshgrid(0.4,linspace(yc, yd, ny));
R_slice3_51 = sqrt((X_slice3_51 - X_slice3_51(1)).^2+(Y_slice3_51-Y_slice3_51(1)).^2);
U_slice3_51 = interp2(x_rect,y_rect,reshape(p,nx,ny)',X_slice3_51,Y_slice3_51,'linear');
save(filename,'X_slice3_51','Y_slice3_51','R_slice3_51','U_slice3_51','-append') 
figure;plot(R_slice3_51,U_slice3_51,'LineWidth',1.5);

% x=0.2
[X_slice4_51, Y_slice4_51] = meshgrid(0.2,linspace(yc, yd, ny));
R_slice4_51 = sqrt((X_slice4_51 - X_slice4_51(1)).^2+(Y_slice4_51-Y_slice4_51(1)).^2);
U_slice4_51 = interp2(x_rect,y_rect,reshape(p,nx,ny)',X_slice4_51,Y_slice4_51,'linear');
save(filename,'X_slice4_51','Y_slice4_51','R_slice4_51','U_slice4_51','-append') 
figure;plot(R_slice4_51,U_slice4_51,'LineWidth',1.5);


