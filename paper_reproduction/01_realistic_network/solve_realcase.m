function solve_realcase(nx,ny)
root=fileparts(mfilename('fullpath'));
addpath(fileparts(root));setup_paths;
xa=0;xb=700;yc=0;yd=600;

total_clock=tic;

x=linspace(xa,xb,nx+1); x=x(:); sx = numel(x)-1; %sx为x方向单元剖分数
y=linspace(yc,yd,ny+1); y=y(:); sy = numel(y)-1; %sy为y方向单元剖分数
celldim = [sx,sy]; %获取网格单元维度
numC = sx*sy;  %单元数
numN = (sx+1)*(sy+1); %节点数
numFX = (sx+1)*sy; numFY = sx*(sy+1); %y方向上面数量； x方向上面数量
numF = numFX + numFY; %面数

[xCoords, yCoords] = ndgrid(x,y);
coords = [xCoords(:), yCoords(:)]; %节点坐标，按照x方向进行编号

N = reshape(1:numN,[sx+1, sy+1]);
NF1 = reshape(N(1:sx+1,1:sy),1,[]);
NF2 = reshape(N(1:sx+1,2:sy+1),1,[]);
faceNodesX = reshape([NF1;NF2],[],1); %x方向上的面编号
NF1 = reshape(N(1:sx,1:sy+1),1,[]);
NF2 = reshape(N(2:sx+1,1:sy+1),1,[]);
faceNodesY = reshape([NF2;NF1],[],1);
faceNodes = [faceNodesX;faceNodesY]; %面对应节点编号

foffset = 0;
FX = reshape(foffset+(1:numFX),sx+1,sy);
foffset = foffset + numFX;
FY = reshape(foffset+(1:numFY),sx,sy+1);
F1 = reshape(FX(1:sx,:),1,[]); %West == 4
F2 = reshape(FX(2:sx+1,:),1,[]); %East == 2
F3 = reshape(FY(:,1:sy),1,[]); %South == 1
F4 = reshape(FY(:,2:sy+1),1,[]); %North == 3
cellFaces = [reshape([F3;F2;F4;F1],[],1),kron(ones([numC,1]),[1,2,3,4]')];

C = zeros([sx+2,sy+2]);
C(2:sx+1,2:sy+1) = reshape(1:numC,[sx,sy]);
NX1 = reshape(C(1:sx+1,2:sy+1),[],1); %左临单元编号
NX2 = reshape(C(2:sx+2,2:sy+1),[],1); %右临单元编号
NY1 = reshape(C(2:sx+1,1:sy+1),[],1); %下临单元编号
NY2 = reshape(C(2:sx+1,2:sy+2),[],1); %上临单元编号
neighbors = [[NX1,NX2];[NY1,NY2]];

NX3 = reshape(C(1:sx,2:sy+1),[],1); %左临单元编号
NX4 = reshape(C(3:sx+2,2:sy+1),[],1); %右临单元编号
NY3 = reshape(C(2:sx+1,1:sy),[],1); %下临单元编号
NY4 = reshape(C(2:sx+1,3:sy+2),[],1); %上临单元编号
cell2cell = [NY3,NX4,NY4,NX3]; %南-东-北-西
G.cell2cell = cell2cell;

k=reshape(N(1:nx,1:ny),[],1);
cNodes =[k k+1 k+(sx+2) k+(sx+1)]; %单元四个节点编号，与Mrst中不太一致

G.cells = struct('num',numC,'facePos',(1:4:(numC+1)*4)','indexMap',(1:numC)');
G.faces = struct('num',numF,'nodePos',(1:2:(numF+1)*2)','neighbors',neighbors,'tag',zeros(numF,1));
G.nodes = struct('num',numN,'coords',coords);
G.cells.faces = cellFaces; G.faces.nodes = faceNodes; G.cellNodes = cNodes;
G.cartDims = celldim; G.griddim = numel(G.cartDims);

edges = reshape(G.faces.nodes,2,[]).'; %获取每条边两个节点编号
[n1, n2] = deal(G.nodes.coords(edges(:,1),:), G.nodes.coords(edges(:,2),:)); %获取对应的坐标
edgeLength = n2 - n1; %获取每条边向量
faceAreas = sqrt(sum(edgeLength.^2,2)); %每条边面积
faceCentroids = (n1+n2)./2; %每条边中心
faceNormals = [edgeLength(:,2), -edgeLength(:,1)]; %每条边法向量，长度为边长

numfaces = diff(G.cells.facePos); %计算每个单元包含边的数量
c = faceCentroids(G.cells.faces(:,1),:); %获取每个单元四条面的面中心坐标
w=1; n=numfaces; %参数
no = rldecode(1:numel(n),n,2).'; %获取每条边所属单元编号
accum = sparse(no,1:numel(no),w); c = accum*[c,ones([size(c,1),1])]; %accum构造每个单元与面对应矩阵，c即表示四个面中心坐标之和
w = c(:,end); c= bsxfun(@rdivide,c(:,1:end-1),w);%w为权重。c为四个面中心坐标之和除以w得到单元中心坐标
cCenter = c; cellno = no; %单元中心坐标。每条边属于哪个单元

cellEdges = edges(G.cells.faces(:,1),:); %获取一个单元四条边每条边的编号大小为[4*numC,2]
r = G.faces.neighbors(G.cells.faces(:,1),2) == cellno; %修正两个节点编号
cellEdges(r,:) = cellEdges(r,[2,1]);
cc = cCenter(cellno,:);
a = G.nodes.coords(cellEdges(:,1),:) - cc;
b = G.nodes.coords(cellEdges(:,2),:) - cc;
quadArea = @(a,b)abs(a(:,1).*b(:,2)-a(:,2).*b(:,1)); %面积计算公式
subArea = quadArea(a,b)./2; %子单元面积
subCentroid = (cCenter(cellno,:) + 2*faceCentroids(G.cells.faces(:,1),:))/3; %子单元中心

c = subCentroid; w=subArea; n = numfaces;
no = rldecode(1:numel(n),n,2).';
accum = sparse(no,1:numel(no),w); c = accum*[c,ones([size(c,1),1])];
w = c(:,end); c= bsxfun(@rdivide,c(:,1:end-1),w);
cellCentroids = c; cellVolumes = w;

G.faces.areas = faceAreas;
G.faces.normals = faceNormals;
G.faces.centroids = faceCentroids;

G.cells.volumes = cellVolumes;
G.cells.centroids = cellCentroids;

zface = G.faces.centroids(:,1); %即查找边界中点x坐标
maxz = max(zface); minz = min(zface);
rightfaces = find(G.faces.centroids(:,1) > (maxz - 1e-6));
leftfaces = find(G.faces.centroids(:,1) < (minz + 1e-6));

nrf = numel(rightfaces); %右侧面Dirichlet
nlf = numel(leftfaces);  %左侧面Dirichlet
bc = struct('face',[],'type',{{}},'value',[]);
bc.face = [leftfaces; rightfaces]; %边界面编号
bc.type = repmat({'pressure'},[1,numel(bc.face)]);
value = zeros(nrf+nlf,1); value(1:nlf) = 1013250; %左边界值1013250 pa，右边界值0
bc.value =value;

[nc,nf] = deal(G.cells.num, G.faces.num);
state = struct('pressure',zeros([nc,1]), 'flux',zeros([nf,1])); %初始化压力与通量解

perm = 1e-14;
perm = repmat(perm,G.cells.num,1);
rock = struct('perm',perm); %渗透率常数

cellno = rldecode(1:G.cells.num,diff(G.cells.facePos),2).'; %获取每个面对应的单元编号
cellNo = cellno;
C = G.cells.centroids; %获取每个单元中心坐标
C = G.faces.centroids(G.cells.faces(:,1),:) - C(cellNo,:); %每个单元每条边中心坐标减去但单元中心坐标

cf = G.cells.faces(:,1); %获取每个单元四条边的编号
sgn = 2*(cellNo == G.faces.neighbors(cf,1)) - 1;  %判断单元四条边是否为流入或者流出
N = bsxfun(@times,sgn,G.faces.normals(cf,:)); %计算外法向量

K = rock.perm * [1, 0, 1]; K = K(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; %张量计算的位置信息

T = zeros(size(cellNo));
for k=1:size(i,2)
    T = T + (C(:,i(k)).*K(cellNo,k).*N(:,j(k)));
end
T = T./sum(C.*C,2);

is_neg = T<0;
if any(is_neg)
    T(is_neg) = -T(is_neg);
end

tic

FractureData = ...
[1,269.611206,152.05243,356.9240112,310.14123;
2,249.5117187,514.990780001,272.218872,470.97082;
3,258.3590698,515.574580001,271.9851684,490.9682;
4,270.6622924,524.702640001,269.1347046,147.78143;
5,355.8302002,348.479800001,337.5810733205,600;
6,366.9730835,338.132990001,426.9185141723,600;
7,198.237915,222.724420001,175.1561889,597.603030001;
8,151.2785034,261.724610001,154.4623059774,600;
9,29.5026855,300.724610001,96.3599853,514.82739;
10,386.0808105,33.3621800002,440.585083,275.191830001;
11,459.6350708,40.2413900001,461.751709,204.812620001;
12,297.180603,237.62103,468.1018066,40.2413900001;
13,312.5264892,272.01678,417.3016967,140.7832;
14,330.5181884,298.47522,439.5266723,156.6582;
15,340.5723877,320.70019,367.5598755,286.304380001;
16,492.9725952,312.762820001,576.5811157,419.6546;
17,505.6726684,309.05859,576.0520019,405.367190001;
18,537.4227905,297.94598,623.3187866,376.68463;
19,322.5338745,380.76941,521.8778076,593.552180001;
20,344.9320678,481.56122,409.8867798,503.959410001;
21,371.8098755,468.12219,510.6787109,383.009210001;
22,432.2849731,510.678830001,642.8280029,374.04999;
23,527.528634971,600,700,473.015615092;
24,0,333.73321,441.2443847,0;
25,13.4389038,342.692380001,347.171875,595.791990001;
26,22.3981933,450.203790001,311.3347778,291.176630001;
27,26.8778076,506.199220001,199.343811,400.92779;
28,44.7963867,528.597410001,365.0905151,342.692380001;
29,378.5294189,309.095210001,512.918518,116.470640001;
30,461.4027099,253.099610001,530.8370971,134.38922;
31,347.171875,374.04999,640.5881958,253.099610001;
32,490.5203857,268.77844,564.4343872,145.58844;
33,47.0361938,181.425410001,53.7556152,306.85541;
34,382.4152832,424.151000001,447.8997192,371.76343;
35,587.9967651,394.78222,549.1029663,362.635190001;
36,589.9812011,393.59161,527.6716919,313.8194;
37,597.125,378.90722,533.6248169,295.960200001;
38,533.6248169,448.75738,453.8527832,326.91638;
39,511.7966919,461.85419,489.5715942,395.17901;
40,565.3748779,425.34161,483.6184692,315.40698;
41,534.4185791,407.482240001,467.3466186,315.803830001;
42,627.2874756,527.3388,574.8999023,498.763610001;
43,644.3532104,519.00439,586.4093017,490.03241;
44,655.8626098,502.335630001,602.6812133,476.53863;
45,415.355896,585.679380001,391.9401855,561.47003;
46,417.3402099,578.535580001,397.8933105,554.326230001;
47,403.0526733,592.029420001,382.0183105,561.86682;
48,495.1278686,505.113580001,468.1403198,481.30121;
49,533.6248169,254.84381,420.9121093,159.196590001;
50,508.6217041,221.10943,441.152771,159.59363;
51,418.5308838,229.04681,312.961914,93.3154300004;
52,362.5714111,174.6748,322.883789,120.69983;
53,357.8088989,216.3468,295.102478,114.74658;
54,402.2589111,283.41882,366.1433105,226.66559;
55,337.5681762,253.256220001,374.4776001,211.18744;
56,386.7808838,264.765620001,509.8123169,101.25281;
57,473.2996826,278.65643,561.0092163,144.909240001;
58,471.7122192,253.653200001,554.6593017,129.034240001;
59,559.0249023,219.125,567.3593139,153.64044;
60,567.7561035,214.759400001,573.7092895,162.37182;
61,574.8999023,215.553040001,579.6624145,173.88104;
62,557.0404663,285.006410001,600.6968994,325.48761;
63,565.3748779,283.022030001,607.0468139,323.503230001];
FractureData(:,1)=[]; FractureData = (FractureData-300)*(1-1e-10)+300;

fractures=zeros(size(FractureData,1),11);
fractures(:,1)=(FractureData(:,1)+FractureData(:,3))/2;
fractures(:,2)=(FractureData(:,2)+FractureData(:,4))/2;
fractures(:,3)=atan2(FractureData(:,4)-FractureData(:,2),FractureData(:,3)-FractureData(:,1));
fractures(:,4)=sqrt((FractureData(:,1)-FractureData(:,3)).^2+(FractureData(:,2)-FractureData(:,4)).^2);
fractures(:,5)=1e-2; %厚度
fractures(:,6)=1e-8; %渗透率
fractures(:,7:10)=[FractureData(:,1),FractureData(:,3),FractureData(:,2),FractureData(:,4)]; %xa,xb,ya,yb

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

fractures(:,11)=0;  %经过的单元数
fractures(:,12)=0;  %经过的总面数
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1); %裂缝开始单元
    preciousfracelem=currentfracselem;  %初始化裂缝前一单元
    while(currentfracselem ~= startendfracs(i,2)) %裂缝终止单元
        [~,~,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1); %相交边界编号
        neighborfracsindex = Cell2cell(currentfracselem,edgenum); %裂缝所在下一个单元,当裂缝与多个边界相交时，需要排除一个重复的编号
        nextfracselem = setdiff(neighborfracsindex,preciousfracelem); %按照单元编号排除，并非按照局部边编号排除
        preciousfracelem = currentfracselem;
        currentfracselem = nextfracselem;
        fractures(i,11)=fractures(i,11)+1;
        fractures(i,12)=fractures(i,12)+size(edgenum,1);
    end

    [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1);
    fractures(i,11)=fractures(i,11)+1;
    fractures(i,12)=fractures(i,12)+size(edgenum,1);
end

pathfractures = zeros(sum(fractures(:,11)),1);
lengthfractures = zeros(sum(fractures(:,11)),1);
centerfractures = zeros(sum(fractures(:,11)),2);   %裂缝中心坐标
currentelemindex = 0;
v_vel = [];
len_fracs = [];
for i=1:size(fractures,1)
    currentfracselem = startendfracs(i,1);
    currentelemindex = currentelemindex + 1;
    pathfractures(currentelemindex) = currentfracselem;
    previousedgenum=0;
    v=[cos(fractures(i,3)), sin(fractures(i,3))];
    v_vel = [v_vel;[repmat(fractures(i,5)*fractures(i,6),fractures(i,11).*size(G.cellNodes,2),1),repmat([cos(fractures(i,3)), sin(fractures(i,3))],fractures(i,11).*size(G.cellNodes,2),1)]];
    while(currentfracselem ~= startendfracs(i,2)) %非裂缝所在结束单元
        [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
            fractures(i,[7,8]),fractures(i,[9,10]));
        edgenum = ii(:,1); %相交边界编号
        if (size(edgenum,1) == 1)
            lengthfractures(currentelemindex) = sqrt((xi-fractures(i,7))^2+(yi-fractures(i,9))^2);
        elseif (size(edgenum,1) == 2)
            lengthfractures(currentelemindex) = sqrt(diff(xi)^2+diff(yi)^2);
        end
        tmp_len = lengthfractures(currentelemindex);
        len_fracs = [len_fracs;[repmat(tmp_len,size(G.cellNodes,2),1),repmat([sin(fractures(i,3));cos(fractures(i,3))],size(G.cellNodes,2)/2,1)]];
        neighborfracsindex=setdiff(edgenum,previousedgenum);
        nextfracselem = Cell2cell(currentfracselem,neighborfracsindex); %裂缝所在下一个单元,当裂缝与多个边界相交时，需要排除一个重复的编号
        previousedgenum = find(Cell2cell(nextfracselem,:) == currentfracselem); %前一个单元边界索引
        currentfracselem = nextfracselem; %迭代现在单元
        leftfracs = ([xi(neighborfracsindex == edgenum);yi(neighborfracsindex == edgenum)]); %现单元左，下位置信息
        currentelemindex = currentelemindex + 1;  %记录此时的
        pathfractures(currentelemindex) = currentfracselem;
    end

    [xi,yi,ii] = polyxpoly(Nodes(Elems(currentfracselem,[1;2;3;4;1]),1),Nodes(Elems(currentfracselem,[1;2;3;4;1]),2),...
        fractures(i,[7,8]),fractures(i,[9,10]));
    edgenum = ii(:,1);
    lengthfractures(currentelemindex) = sqrt((xi-fractures(i,8))^2+(yi-fractures(i,10))^2);
    tmp_len = lengthfractures(currentelemindex);
    len_fracs = [len_fracs;[repmat(tmp_len,size(G.cellNodes,2),1),repmat([sin(fractures(i,3));cos(fractures(i,3))],size(G.cellNodes,2)/2,1)]];

end

time2 = toc
posfracs1 = (pathfractures(:,1)'-1)*size(G.cellNodes,2)+[1:4]'; posfracs1 = posfracs1(:);
facsfracs1 = G.cells.faces(posfracs1(:),1);
uN_fracs2 = N(posfracs1,:)./sqrt(N(posfracs1,1).^2+N(posfracs1,2).^2);

posfracs22 = (pathfractures(:,1)'-1)*size(G.cellNodes,2)+[2 3 4 1]';
posfracs22 = posfracs22(:);
facsfracs2 = G.cells.faces(posfracs22,1);

v_p_n = v_vel(:,2).*uN_fracs2(:,1)+v_vel(:,3).*uN_fracs2(:,2);
active_edges = find(v_p_n ~= 0);
posfracs1 = posfracs1(active_edges,:);
facsfracs1 = facsfracs1(active_edges,:);
v_vel = v_vel(active_edges,:);
uN_fracs2 = uN_fracs2(active_edges,:);

kfvn = 1./abs(v_vel(:,2).*uN_fracs2(:,1)+v_vel(:,3).*uN_fracs2(:,2));
kfvvt = kfvn.*v_vel(:,1).* ...
    [v_vel(:,2).*v_vel(:,2),v_vel(:,2).*v_vel(:,3),v_vel(:,3).*v_vel(:,3)];
K_Fracs = kfvvt(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; %张量计算的位置信息
C_fracs = C(posfracs1,:);
N_fracs = uN_fracs2;

w_fracs =len_fracs(:,1).*len_fracs(:,2)./G.faces.areas(facsfracs2,:);
w_fracs = abs(w_fracs);
C_fracs = C_fracs.*w_fracs;

T_fracs = zeros(size(facsfracs1,1),1);
for k=1:size(i,2)
    T_fracs = T_fracs + (C_fracs(:,i(k)).*K_Fracs(:,k).*N_fracs(:,j(k)));
end
T_fracs = T_fracs./sum(C_fracs.*C_fracs,2);

posfracs = posfracs1;

T_matrix = T; % diagnostic copy before adding fracture transmissibilities
T = T + sparse(posfracs,ones(size(posfracs)),T_fracs,size(T,1),size(T,2));

neighborship = G.faces.neighbors;   nif = size(neighborship,1); %每条边与哪些单元相邻
cf = G.cells.faces(:,1);  ncf = size(cf,1); %每个单元四条边全局边编号
nc = G.cells.num;  n = nc; %全局自由度数量
ft = 1./accumarray(cf,1./T,[nif,1]);  %调和平均计算面传导率系数

if ~isempty(bc)
    hh=zeros(nif,1);  %Neumann边界值
    ff=zeros(ncf,1);  %Dirichlet边界值
    gg=zeros(nc,1);   %源项值

    dF = false([G.faces.num, 1]);
    dC = [];
    is_press = strcmpi('pressure',bc.type);
    face = bc.face(is_press); %获得Dirichlet面全局编号
    dC = bc.value(is_press); %获得Dirichlet面边界值
    map = sparse(double(face),1,1:numel(face)); %为后面每个面赋值做映射
    dF(face) = true;  %将Dirichlet边界条件判断分配给全局面编号中
    i = dF(G.cells.faces(:,1)); %将其分配给单元对应四个面矩阵中
    ff(i) = -dC(map(G.cells.faces(i,1))); %赋值给ff中
    dC = dC(map(dF));  %按照dF顺序重新排列赋值

    is_flux = strcmpi('flux',bc.type);
    hh(bc.face(is_flux)) = -bc.value(is_flux);

end

src = [];
if ~isempty(src)
    ss = accumarray(src.cell.src.value)
    ii = accumarray(src.cell,1)>0;
    gg(ii) = gg(ii) +ss(ii);
end

i=all(neighborship ~=0, 2); %获取内部边界
j = i(cf) | dF(cf); %根据cf即每个单元包含四条边来判断边界是否为非Neumann边界

rhs = accumarray(cellNo,-ft(cf).*ff,[n,1]) + gg  + accumarray(cellNo,-hh(cf),[n,1]);

d = zeros(G.cells.num,1);
d = d + accumarray(cellNo(dF(cf)),T(dF(cf)),[nc,1]) +...
    accumarray(reshape(neighborship(i,:),[],1), repmat(ft(i),[2,1]), [nc,1]);

I = [neighborship(i,1); neighborship(i,2); (1:nc)'];
J = [neighborship(i,2); neighborship(i,1); (1:nc)'];
V = [-ft(i); -ft(i); d];
A = sparse(double(I), double(J), V, nc, nc);
p = A\rhs; %求解全局系统

fpress = accumarray(cf, p(cellNo).*T, [nif,1])./accumarray(cf(:,1),T,[nif,1]);

b = any(G.faces.neighbors==0,2);
fpress(b) = fpress(b) - hh(b)./ft(b);

fpress(dF)=dC;

sgn = 2*(G.faces.neighbors(~i,2)==0)-1;
ni = neighborship(i,:); %内部边单元相邻关系
flux = -accumarray(find(i), ft(i).*(p(ni(:,2))-p(ni(:,1))),[nif,1]); %内部边通量
c = sum(G.faces.neighbors(~i,:),2); %边界边相邻单元
flux(~i) = -sgn.*ft(~i).*(fpress(~i) - p(c)); %由公式得到边界边通量

state.pressure(1:nc) = p(1:nc);
state.flux = flux;
state.facePressure = fpress;

seconds=toc(total_clock);
cellNo=repelem((1:G.cells.num)',diff(G.cells.facePos));
faceNo=G.cells.faces(:,1);
sgn=2*(cellNo==G.faces.neighbors(faceNo,1))-1;
residual=accumarray(cellNo,sgn.*state.flux(faceNo),[G.cells.num,1])-gg;
fractured=false(G.cells.num,1);
fractured(unique(pathfractures))=true;
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'mass_balance_data.mat'),'G','state','fractures','fractured','residual','nx','ny','seconds');
writetable(table((1:G.cells.num)',G.cells.centroids(:,1),G.cells.centroids(:,2),fractured,residual,abs(residual), ...
    'VariableNames',{'Cell','x','y','Fractured','SignedResidual','AbsoluteResidual'}),fullfile(out,'residuals.csv'));
fprintf('Cells: %d; fractured: %d; maximum absolute residual: %.8e\n',G.cells.num,sum(fractured),max(abs(residual)));
end
function A = rldecode(A, n, dim)

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

d      = 1:max(dim, ndims(A));
d([1, dim])   = [dim, 1];
B      = permute(A,d);

r      = n(:)~=0;
B      = reshape(B(r, :), sum(r), []);

i      = cumsum([1; double(reshape(n(r), [], 1))]);
j      = zeros(i(end)-1,1);
j(i(1:end-1)) = 1;

szA    = [size(A), ones(1, dim-ndims(A))];
A      = permute(reshape(B(cumsum(j),:), [sum(n(:)), szA(d(2:end))]), d);
end
