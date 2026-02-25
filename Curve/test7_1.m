%curve fracture test
%在每个单元上基于采样点仅仅存在一个点，化曲为直
%相关边界条件，渗透率，厚度等与Cross case 几乎保持一致。

%圆形裂缝 Circle

clc,clear
close all

nnn=20;
nx=nnn; ny=nnn;
xa=0;xb=1;yc=0;yd=1;
hx=(xb-xa)/nx;
hy=(yd-yc)/ny;
tic

%% 网格信息
x=linspace(xa,xb,nx+1); x=x(:); sx = numel(x)-1; %sx为x方向单元剖分数
y=linspace(yc,yd,ny+1); y=y(:); sy = numel(y)-1; %sy为y方向单元剖分数
celldim = [sx,sy]; %获取网格单元维度  
numC = sx*sy;  %单元数
numN = (sx+1)*(sy+1); %节点数
numFX = (sx+1)*sy; numFY = sx*(sy+1); %y方向上面数量； x方向上面数量
numF = numFX + numFY; %面数

% Nodes/Coordinates 
[xCoords, yCoords] = ndgrid(x,y);
coords = [xCoords(:), yCoords(:)]; %节点坐标，按照x方向进行编号

% Generate face-edges
N = reshape(1:numN,[sx+1, sy+1]);
NF1 = reshape(N(1:sx+1,1:sy),1,[]); 
NF2 = reshape(N(1:sx+1,2:sy+1),1,[]);
faceNodesX = reshape([NF1;NF2],[],1); %x方向上的面编号
NF1 = reshape(N(1:sx,1:sy+1),1,[]);
NF2 = reshape(N(2:sx+1,1:sy+1),1,[]);
faceNodesY = reshape([NF2;NF1],[],1);
faceNodes = [faceNodesX;faceNodesY]; %面对应节点编号

% Generate cell-faces
foffset = 0;
FX = reshape(foffset+(1:numFX),sx+1,sy);
foffset = foffset + numFX;
FY = reshape(foffset+(1:numFY),sx,sy+1);
F1 = reshape(FX(1:sx,:),1,[]); %West == 4
F2 = reshape(FX(2:sx+1,:),1,[]); %East == 2
F3 = reshape(FY(:,1:sy),1,[]); %South == 1
F4 = reshape(FY(:,2:sy+1),1,[]); %North == 3
%每个单元对应面编号，第一列为面编号，第二列为局部编号
cellFaces = [reshape([F3;F2;F4;F1],[],1),kron(ones([numC,1]),[1,2,3,4]')]; 

% Generate neighbors
C = zeros([sx+2,sy+2]);
C(2:sx+1,2:sy+1) = reshape(1:numC,[sx,sy]);
NX1 = reshape(C(1:sx+1,2:sy+1),[],1); %左临单元编号
NX2 = reshape(C(2:sx+2,2:sy+1),[],1); %右临单元编号
NY1 = reshape(C(2:sx+1,1:sy+1),[],1); %下临单元编号
NY2 = reshape(C(2:sx+1,2:sy+2),[],1); %上临单元编号
neighbors = [[NX1,NX2];[NY1,NY2]];

%新增 cells neighbors cells
NX3 = reshape(C(1:sx,2:sy+1),[],1); %左临单元编号
NX4 = reshape(C(3:sx+2,2:sy+1),[],1); %右临单元编号
NY3 = reshape(C(2:sx+1,1:sy),[],1); %下临单元编号
NY4 = reshape(C(2:sx+1,3:sy+2),[],1); %上临单元编号
cell2cell = [NY3,NX4,NY4,NX3]; %南-东-北-西
G.cell2cell = cell2cell;

% Generate cell nodes
k=reshape(N(1:nx,1:ny),[],1);
cNodes =[k k+1 k+(sx+2) k+(sx+1)]; %单元四个节点编号，与Mrst中不太一致

% Assemble structure
G.cells = struct('num',numC,'facePos',(1:4:(numC+1)*4)','indexMap',(1:numC)');
G.faces = struct('num',numF,'nodePos',(1:2:(numF+1)*2)','neighbors',neighbors,'tag',zeros(numF,1));
G.nodes = struct('num',numN,'coords',coords);
G.cells.faces = cellFaces; G.faces.nodes = faceNodes; G.cellNodes = cNodes;
G.cartDims = celldim; G.griddim = numel(G.cartDims);

%补充相关信息
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

% Update grid
G.faces.areas = faceAreas;
G.faces.normals = faceNormals;
G.faces.centroids = faceCentroids;

G.cells.volumes = cellVolumes;
G.cells.centroids = cellCentroids;

%% 初始化渗透率、边界条件与解
zface = G.faces.centroids(:,1); %即查找边界中点x坐标
maxz = max(zface); minz = min(zface);
rightfaces = find(G.faces.centroids(:,1) > (maxz - 1e-6));
leftfaces = find(G.faces.centroids(:,1) < (minz + 1e-6));

nrf = numel(rightfaces); %右侧面Dirichlet
nlf = numel(leftfaces);  %左侧面Dirichlet
bc = struct('face',[],'type',{{}},'value',[]);
bc.face = [leftfaces; rightfaces]; %边界面编号
bc.type = repmat({'pressure'},[1,numel(bc.face)]);
value = zeros(nrf+nlf,1); value(1:nlf) = 1; %左边界值1，右边界值0
bc.value =value;

[nc,nf] = deal(G.cells.num, G.faces.num);
state = struct('pressure',zeros([nc,1]), 'flux',zeros([nf,1])); %初始化压力与通量解

%基质单元每个单元渗透率为1
perm = 1;
perm = repmat(perm,G.cells.num,1);
rock = struct('perm',perm); %渗透率常数

%% 计算 T
%这一部分是计算传导率系数
cellno = rldecode(1:G.cells.num,diff(G.cells.facePos),2).'; %获取每个面对应的单元编号
cellNo = cellno;
C = G.cells.centroids; %获取每个单元中心坐标
C = G.faces.centroids(G.cells.faces(:,1),:) - C(cellNo,:); %每个单元每条边中心坐标减去但单元中心坐标

cf = G.cells.faces(:,1); %获取每个单元四条边的编号
sgn = 2*(cellNo == G.faces.neighbors(cf,1)) - 1;  %判断单元四条边是否为流入或者流出
N = bsxfun(@times,sgn,G.faces.normals(cf,:)); %计算外法向量

%将K扩展为张量
K = rock.perm * [1, 0, 1]; K = K(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; %张量计算的位置信息

%计算传导率系数 % Compute T = C'*K*N / C'*C.
T = zeros(size(cellNo));
for k=1:size(i,2)
    T = T + (C(:,i(k)).*K(cellNo,k).*N(:,j(k)));
end
T = T./sum(C.*C,2);

%修正负传导率系数
is_neg = T<0;
if any(is_neg)
    T(is_neg) = -T(is_neg);
end

%% RDFM裂缝 求解T test2 对四条边均有影响
%裂缝相关信息，行代表裂缝数量，列代表裂缝具体信息
%圆形裂缝
numfractures=0;
numfractures=numfractures+1;
% list of large prime numbers used to choose N :
% 1009,2003,3001,5003,10007,20011,30011,40009,50021,60013,70001,80021,90001,100003
N_partition = 100003; 
t_start = 0; t_end = 2*pi;
t_partition = linspace(t_start,t_end,N_partition)';
Radius = 1/3;%半径
x_center = 0.5+(rand-0.5)*1e-4; y_center = 0.5+(rand-0.5)*1e-4;  %裂缝参数，（x_center,y_center）即为中心坐标
x_partition = x_center + Radius*cos(t_partition); y_partition = y_center + Radius*sin(t_partition); %裂缝参数方程
%判断裂缝位于哪个单元内部，由于单元采用横向编号，故我们采用以下方法可判定出裂缝大致位于哪个单元内部中
figure;
%节点坐标
x=linspace(xa,xb,nx+1);
y=linspace(yc,yd,ny+1);
[xx,yy]=meshgrid(x,y); %横向编号
surf(xx,yy,zeros(nx+1,ny+1));
hold on; plot(x_partition,y_partition,'r')
ele_located = (ceil((y_partition-yc)/hy)-1)*nx + ceil((x_partition-xa)/hx);
where_inter=find(diff([-1;ele_located])~=0);   %判断在哪个单元，只记录首次出现，负一是记录首个单元
ele_passed = ele_located(where_inter);
numberoffracturespassed(numfractures,1) = length(ele_passed);
pathfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + [1:length(ele_passed)],1)=ele_passed;  %裂缝经过单元,为包含多组弯曲裂缝做设计

%新增，也不算新增，离散化表示积分点上裂缝相关信息
%我们省略高斯点，直接带入端点进行模拟
permeabilityfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + [1:length(ele_passed)],2) = 0; %裂缝渗透率系数
thicknessfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + [1:length(ele_passed)],2) = 0; %裂缝厚度
lengthfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + [1:length(ele_passed)],1) = 0;

Nodes = G.nodes.coords;
Elems = G.cellNodes';
Cell2cell = G.cell2cell';

%计算每个单元与裂缝的两个交点
frac_nodes=zeros(2*length(ele_passed),2); %第一列表示x坐标，第二列表示y坐标， 第[2*n-1,2*n]行表示与单元的两个交点
theta_avg_sum = zeros(length(ele_passed),1);
v = zeros(length(ele_passed),2); %第一列表示cos(theta_avg)，第二列表示sin(theta_avg)；
for ell = 1:length(ele_passed)
    %判断裂缝经过的边
    if ell==1  %弯曲裂缝经过初始单元
        edge_in  = -1;
        edge_out = find(Cell2cell(:,ele_passed(ell)) == ele_passed(ell+1));
        time_in = t_start;
    elseif ell == length(ele_passed) %弯曲裂缝经过终止单元
        edge_in = find(Cell2cell(:,ele_passed(ell)) == ele_passed(ell-1));
        edge_out = -1;
        time_out =t_end;
    else   %其余裂缝经过单元
        edge_in = find(Cell2cell(:,ele_passed(ell)) == ele_passed(ell-1));
        edge_out = find(Cell2cell(:,ele_passed(ell)) == ele_passed(ell+1));
    end
    
    if(edge_in == 1) || (edge_in == 3) %1即为bottom边界，3即为top边界
        theta = asin(((Nodes(Elems(1,ele_passed(ell)),2)+(edge_in-1)/2*hy)-y_center)/Radius); %计算角度
        t_candidates = [theta;pi-theta;2*pi+theta]; %???
        time_in = t_candidates(t_candidates>=t_partition(where_inter(ell)-1) & t_candidates<= t_partition(where_inter(ell)));  %计算进入此单元的t参数信息
    end
    
    if(edge_in == 2) || (edge_in == 4) %2即为right边界，4即为left边界
        theta = acos(((Nodes(Elems(1,ele_passed(ell)),1)+(4-edge_in)/2*hx)-x_center)/Radius); %计算角度
        t_candidates = [theta;2*pi-theta]; %???
        time_in = t_candidates(t_candidates>=t_partition(where_inter(ell)-1) & t_candidates<= t_partition(where_inter(ell)));  %计算进入此单元的t参数信息
    end
    
    if (edge_out == 1) || (edge_out == 3)
        theta = asin(((Nodes(Elems(1,ele_passed(ell)),2)+(edge_out-1)/2*hy) - y_center) / Radius );
        t_candidates = [theta;pi-theta;2*pi+theta];
        time_out = t_candidates(t_candidates>=t_partition(where_inter(ell+1)-1) & t_candidates<=t_partition(where_inter(ell+1)));
    end
    
    if (edge_out == 2) || (edge_out == 4)
        theta = acos( ((Nodes(Elems(1,ele_passed(ell)),1)+(4-edge_out)/2*hx) - x_center) / Radius );
        t_candidates = [theta;2*pi-theta];
        time_out = t_candidates(t_candidates>=t_partition(where_inter(ell+1)-1) & t_candidates<=t_partition(where_inter(ell+1)));
    end  
    
    %有限体积中，暂时只考虑两个采样点
    t_sampling = time_in + (time_out - time_in)*[0;1];
    x_sampling = x_center + Radius*cos(t_sampling);          %每个单元采样点x坐标
    y_sampling = y_center + Radius*sin(t_sampling);        %每个单元采用点y坐标
    
    frac_nodes(2*ell-1:2*ell,:) = [x_sampling,y_sampling];
    
    %取两点斜率作为theta_avg
    theta_avg1 = atan2((y_sampling(1)-y_sampling(2)),(x_sampling(1)-x_sampling(2)));
    if theta_avg1 < 0
        theta_avg = theta_avg1 + pi;  %修正将转化为[0,pi]内
    else
        theta_avg = theta_avg1;
    end
    
    %设置theta
    theta_avg_sum(ell,:) = theta_avg;
    %设置v
    v(ell,:) = [cos(theta_avg),sin(theta_avg)];

    
    permeabilityfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + ell,1:2) = 1e8 * ones(size(t_sampling));
    thicknessfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + ell,1:2) = 0.005 * ones(size(t_sampling));
    
    lengthfractures(sum(numberoffracturespassed(1:numfractures-1,1)) + ell,1) = time_out - time_in;
    
end


%已知信息，裂缝经过单元，
%裂缝与单元的两个交点，裂缝在每个单元上的长度，两个交点构成直线的倾角
%以下为计算裂缝通量
%虚拟裂缝部分，即增加权重,这一部分与test4_5部分一致

%求解裂缝中点位置坐标，借助得到的frac_nodes进行求解,其中fracs_nodes的行数必为pathfractures的两倍
frac_mid_nodes = squeeze(sum(reshape(frac_nodes, 2, size(frac_nodes,1)/2, size(frac_nodes,2)),1))/2;

%加入随机扰动,用来排除  单元中心的影响
frac_mid_nodes = frac_mid_nodes + (rand(size(frac_mid_nodes,1),size(frac_mid_nodes,2)))*1e-10;
%遍历每个裂缝单元，判定相近边（两条，一条竖线，一条横线）
%定义一个fracs_local_num表示裂缝相近的局部点，暂时不考虑theta，第一列表示横靠近哪条边，第二列表示竖靠近哪条边
fracs_local_num = zeros(size(frac_mid_nodes,1),2);
fracs_local_dist = zeros(size(frac_mid_nodes,1),2);
Elems1 = G.cellNodes; %特殊设置
for ell = 1:length(ele_passed)
    currentelem = pathfractures(ell);
    local_nodes = Nodes(Elems1(currentelem,[1;2;3;4]),:);
    x1 = local_nodes(1,1); y1 = local_nodes(1,2);
    x2 = local_nodes(2,1); y2 = local_nodes(2,2);
    x3 = local_nodes(3,1); y3 = local_nodes(3,2);
    x4 = local_nodes(4,1); y4 = local_nodes(4,2);
    %单元内裂缝中心坐标
    local_frac_nodes = frac_mid_nodes(ell,:);
    x0 = local_frac_nodes(1,1); y0 = local_frac_nodes(1,2);
    %由点至直线的距离公式计算得到裂缝中心点至每个边界的距离
    %求解裂缝中点至局部边1的距离
    distance(1,:) = abs((y1-y2)*x0-(x1-x2)*y0+x1*y2-x2*y1)/sqrt((y1-y2)^2+(x1-x2)^2);
    %求解裂缝中点至局部边3的距离
    distance(3,:) = abs((y3-y4)*x0-(x3-x4)*y0+x3*y4-x4*y3)/sqrt((y3-y4)^2+(x3-x4)^2);
    
    %排除等于角度为pi/2的情况，
    %考虑到浮点误差，我们设置此值为sin(pi/2)。
    if ( abs(sin(theta_avg_sum(ell))) < abs(sin(pi/2)) )
        local1 = find(distance == min(distance([1 3])));
        dist1 = min(distance([1 3]));
    elseif ( abs(sin(theta_avg_sum(ell))) >=  abs(sin(pi/2)) )
        local1 = -1;
        dist1 = -1;
    end
    
    %求解裂缝中点至局部边2的距离
    distance(2,:) = abs((y2-y3)*x0-(x2-x3)*y0+x2*y3-x3*y2)/sqrt((y2-y3)^2+(x2-x3)^2);
    %求解裂缝中点至局部边4的距离
    distance(4,:) = abs((y4-y1)*x0-(x4-x1)*y0+x4*y1-x1*y4)/sqrt((y4-y1)^2+(x4-x1)^2);
    
    %排除角度为0的情况
    %考虑到浮点误差，我们设置此值为cos(0)。
    if ( abs(cos(theta_avg_sum(ell))) < abs(cos(0)) )
        local2 = find(distance == min(distance([2 4])));
        dist2 = min(distance([2 4]));
    elseif ( abs(cos(theta_avg_sum(ell))) >=  abs(cos(0)) )
        local2 = -1;
        dist2 = -1;
    end
    
    %将相近边带入总数组中
    fracs_local_num(ell,:) = [local1 local2];
    fracs_local_dist(ell,:) = [dist1 dist2];
end

%得到局部边以后，在此排除掉能被pi整除的情况,
%在此我们联立裂缝信息，第一列表示裂缝经过单元，第二列表示虚拟裂缝局部边，第三列表示距离，用来计算权重，第四列表示所属裂缝编号
%裂缝局部编号
fracs_local_num1 = reshape(fracs_local_num',[],1);
%裂缝局部距离
fracs_local_dist1 = reshape(fracs_local_dist',[],1);
%裂缝所属单元
pathfractures1 = rldecode(pathfractures,2,1);
%联立信息
actual_fracs_info  = [pathfractures1,fracs_local_num1,fracs_local_dist1];
%排除标记包含‘-1’的行，即排除一些不存在的虚拟裂缝
mask1 = ((actual_fracs_info(:,2) == -1));
actual_fracs_info1 = actual_fracs_info(~mask1,:);

%以下求解虚拟裂缝的单元信息
%借助上述actual_fracs_info1中第一行与第二行得到的单元信息与局部边信息，即判定得到全局边
%首先得到全局边单元唯一编号
fracs_global_mask = (actual_fracs_info1(:,1)-1)*4 + actual_fracs_info1(:,2);
%其次利用单元唯一编号得到全局边编号
fracs_global_num = G.cells.faces(fracs_global_mask,1);
%获得其相邻边单元
fracs_global_elem = G.faces.neighbors(fracs_global_num,:);
%需要将其进行重新排序
%找出相同位置的元素
[n,m] = size(fracs_global_elem);
same_mask  = (actual_fracs_info1(:,1) == fracs_global_elem);
%列位置索引
[~, sort_idx] = sort(~same_mask, 2);
linear_indices = sub2ind([n, m], repmat((1:n)',1,m), sort_idx);
% 重新排列,由此得到第一列即为我们的实际裂缝的单元，第二列为虚拟裂缝的单元
fracs_global_elem1 = fracs_global_elem(linear_indices);

%除了需要知道虚拟裂缝位于的单元外，还需知道虚拟裂缝的局部边编号，
%因此我们需要反向计算
%即根据单元得到四条边全局边单元唯一编号，再有唯一编号得到边全局编号，再按照局部序号排序为4列
%再由每一行的编号与全局边编号进行一一对比，得到局部边索引
%修正，排除掉虚拟裂缝位于边界单元的值，即排除相邻边为编号0，记为mask2，对应虚拟裂缝也需要更改
mask2 = ((fracs_global_elem1(:,2) == 0));                                                                         %new
virtual_fracs_global_elem1_fix = fracs_global_elem1(~mask2,2);                                                    %new
fracs_global_num_fix = fracs_global_num(~mask2);                                                                  %new
face_test = reshape(G.cells.faces((virtual_fracs_global_elem1_fix-1)*4 + [1:4],1),[],4);                          %new
virtual_local_num1 = arrayfun(@(row)find(face_test(row,:) == fracs_global_num_fix(row),1),1:size(face_test,1))';  %new
% face_test = reshape(G.cells.faces((fracs_global_elem1(:,2)-1)*4 + [1:4],1),[],4);
% virtual_local_num1 = arrayfun(@(row)find(face_test(row,:) == fracs_global_num(row),1),1:size(face_test,1))';

%因此我们得到实际裂缝单元信息与邻近边信息，虚拟单元信息与邻近边信息，即
%第一列表示虚拟裂缝单元信息，第二列表示虚拟裂缝邻近边信息，第三列表示相应的实际裂缝编号
virtual_fracs_info1 = [virtual_fracs_global_elem1_fix, virtual_local_num1];         %new
% virtual_fracs_info1 = [fracs_global_elem1(:,2), virtual_local_num1, actual_fracs_info1(:,4)];


%以下为求解实际裂缝与单元的交点，为计算传导率系数做准备
%首先根据我们的理论，即裂缝由达西定律近似表示，因此我们只要计算一些
%首先计算实际裂缝的实际相交边，由局部边可得到实际相交边，
actual_local_nums = [mod(actual_fracs_info1(:,2)-2,4)+1 mod(actual_fracs_info1(:,2),4)+1];
%计算虚拟裂缝的实际相交边，由局部边可得到实际相交边
virtual_local_nums = [mod(virtual_fracs_info1(:,2)-2,4)+1 mod(virtual_fracs_info1(:,2),4)+1];

%将实际裂缝和虚拟裂缝的信息分别对应至实际相交边上
%首先对应实际裂缝,第一列表示单元，第二列表示相交边
actual_local_nums1 = reshape(actual_local_nums',[],1);
actual_fracs_info2 = [rldecode(actual_fracs_info1(:,1),2,1),actual_local_nums1];
%对应虚拟裂缝,第一列表示单元，第二列表示相交边
virtual_local_nums1 = reshape(virtual_local_nums',[],1);
virtual_fracs_info2 = [rldecode(virtual_fracs_info1(:,1),2,1),virtual_local_nums1];

%由前面已经获得裂缝的裂缝信息进行组装，即第一列表示裂缝在该单元的长度，
%第二列表示相应cos、sin值，第三列表示eplison*k_f，第四和五列表示[cos sin]值
actual_fracs_flux_info = [rldecode(lengthfractures(:,1),2,1), reshape(v',[],1),...
    rldecode(thicknessfractures(:,1).*permeabilityfractures(:,1),2,1), rldecode(v,2,1)];
%借助前面mask1，排除不需要的行，即得到
actual_fracs_flux_info1 = actual_fracs_flux_info(~mask1,:);
%第一列表示单元，第二列表示局部边编号，
%第三列表示裂缝在该单元上的长度，第四列表示对应的角度和cos、sin值
%第五列表示渗透率乘以厚度，第六列表示costheta，第七列表示sintheta，
%第八列表示中点与临近边距离，第九列表示所属裂缝编号
actual_fracs_flux_info2 = [actual_fracs_info2, rldecode(actual_fracs_flux_info1,2,1), rldecode(actual_fracs_info1(:,[3]),2,1)];

%以下为求解实际裂缝通量
%实际裂缝全局唯一编号
actual_fracs_pos = (actual_fracs_flux_info2(:,1)-1)*4 + actual_fracs_flux_info2(:,2);
%计算单位外法向量
actual_fracs_N = N(actual_fracs_pos,:)./sqrt(N(actual_fracs_pos,1).^2+N(actual_fracs_pos,2).^2);
%计算1/|v*n|
kfvn = 1./abs(actual_fracs_flux_info2(:,6).*actual_fracs_N(:,1)+actual_fracs_flux_info2(:,7).*actual_fracs_N(:,2));
%计算eplison*k_f*1/|v*n|*v*vt
kfvvt = kfvn.*actual_fracs_flux_info2(:,5).*...
    [actual_fracs_flux_info2(:,6).*actual_fracs_flux_info2(:,6),...
    actual_fracs_flux_info2(:,6).*actual_fracs_flux_info2(:,7),...
    actual_fracs_flux_info2(:,7).*actual_fracs_flux_info2(:,7)];
K_Fracs = kfvvt(:,[1,2,2,3]);
i = [1,1,2,2]; j=[1,2,1,2]; %张量计算的位置信息
%计算单元中心至边中心向量
C_fracs = C(actual_fracs_pos,:);
N_fracs = actual_fracs_N;
%达西定律转化渗透率
actual_fracs_face = G.cells.faces(actual_fracs_pos,1);
w_fracs = actual_fracs_flux_info2(:,3).*actual_fracs_flux_info2(:,4)./G.faces.areas(actual_fracs_face,:);
w_fracs = abs(w_fracs);
C_fracs = C_fracs.*w_fracs;

%计算传导率系数 % Compute T = C'*K*N / C'*C.
T_fracs = zeros(size(actual_fracs_pos,1),1);
for k=1:size(i,2)
    T_fracs = T_fracs + (C_fracs(:,i(k)).*K_Fracs(:,k).*N_fracs(:,j(k)));
end
T_fracs = T_fracs./sum(C_fracs.*C_fracs,2);

%以下为计算虚拟裂缝全局编号
virtual_fracs_pos = (virtual_fracs_info2(:,1)-1)*4+virtual_fracs_info2(:,2);

%以下计算实际裂缝与虚拟裂缝权重
C_fracs_len = sqrt(C(actual_fracs_pos,1).^2+C(actual_fracs_pos,2).^2);
w1 = (C_fracs_len + actual_fracs_flux_info2(:,8))./(2*C_fracs_len);

%虚拟裂缝权重需要重新分配大小
w2 = 1 - w1;
w2 = w2(~rldecode(mask2,2,1),:);                                       %new


% %测试w1,w2
w1=ones(size(w1,1),1);
w2=zeros(size(w2,1),1);

w=[w1;w2];

%求解总裂缝位置和通量
Fracs_pos_all = [actual_fracs_pos;virtual_fracs_pos];
% T_fracs_all = w.*[T_fracs;T_fracs];
T_fracs_all = w.*[T_fracs;T_fracs(~rldecode(mask2,2,1),:)];            %new

%% 叠加总T
T = T + sparse(Fracs_pos_all,ones(size(Fracs_pos_all)),T_fracs_all,size(T,1),size(T,2));


%% 计算全局系统并求解
%以下计算每条边的传导率系数并组装全局矩阵
neighborship = G.faces.neighbors;   nif = size(neighborship,1); %每条边与哪些单元相邻
cf = G.cells.faces(:,1);  ncf = size(cf,1); %每个单元四条边全局边编号
nc = G.cells.num;  n = nc; %全局自由度数量
ft = 1./accumarray(cf,1./T,[nif,1]);  %调和平均计算面传导率系数

if ~isempty(bc)
    %计算边界条件，为组装全局系统做准备
    hh=zeros(nif,1);  %Neumann边界值
    ff=zeros(ncf,1);  %Dirichlet边界值
    gg=zeros(nc,1);   %源项值
    
    %处理Dirichlet边界
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
    
    %处理Neumann边界
    is_flux = strcmpi('flux',bc.type);
    hh(bc.face(is_flux)) = -bc.value(is_flux);
    
end

src = [];
%源项处理
if ~isempty(src)
    ss = accumarray(src.cell.src.value)
    ii = accumarray(src.cell,1)>0;
    gg(ii) = gg(ii) +ss(ii);
end

i=all(neighborship ~=0, 2); %获取内部边界
j = i(cf) | dF(cf); %根据cf即每个单元包含四条边来判断边界是否为非Neumann边界

%右端项，包含Dirichlet边界处理，源项处理，Neumann边界处理
rhs = accumarray(cellNo,-ft(cf).*ff,[n,1]) + gg  + accumarray(cellNo,-hh(cf),[n,1]);

%左侧矩阵，内部边部分传导率以及Dirichlet边压力相对应的传导率
d = zeros(G.cells.num,1);
d = d + accumarray(cellNo(dF(cf)),T(dF(cf)),[nc,1]) +...
    accumarray(reshape(neighborship(i,:),[],1), repmat(ft(i),[2,1]), [nc,1]); %计算每个单元的总共贡献

%组装矩阵
I = [neighborship(i,1); neighborship(i,2); (1:nc)'];
J = [neighborship(i,2); neighborship(i,1); (1:nc)'];
V = [-ft(i); -ft(i); d];
A = sparse(double(I), double(J), V, nc, nc);
p = A\rhs; %求解全局系统

%构造面压力
fpress = accumarray(cf, p(cellNo).*T, [nif,1])./accumarray(cf(:,1),T,[nif,1]);

%处理Neumann边面压力
b = any(G.faces.neighbors==0,2);
fpress(b) = fpress(b) - hh(b)./ft(b);

%处理Dirichlet边面压力
fpress(dF)=dC;

%计算通量
sgn = 2*(G.faces.neighbors(~i,2)==0)-1;
ni = neighborship(i,:); %内部边单元相邻关系
flux = -accumarray(find(i), ft(i).*(p(ni(:,2))-p(ni(:,1))),[nif,1]); %内部边通量
c = sum(G.faces.neighbors(~i,:),2); %边界边相邻单元
flux(~i) = -sgn.*ft(~i).*(fpress(~i) - p(c)); %由公式得到边界边通量

%总结
state.pressure(1:nc) = p(1:nc);
state.flux = flux;
state.facePressure = fpress;

%% postprocessing 绘图

%绘制网格
% figure;
[x_rect, y_rect] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));

% contourf(x_rect,y_rect,reshape(p,nx,ny)',20,'linestyle','no');
% hold on
% plot(fractures(:,[7,8])',fractures(:,[9,10])',...
%    'r-','LineWidth',1.5); hold off;
% colorbar;
% % title('fractures tpfa');

figure;
[x_rect1, y_rect1] = meshgrid(linspace(xa, xb, nx), linspace(yc, yd, ny));
node1=[x_rect1(:),y_rect1(:)];
surf(x_rect,y_rect,reshape(p,nx,ny)');
Interp_xmin = xa+1e-7; Interp_xmax =xb-1e-7 ; Interp_ymin = yc+1e-7; Interp_ymax = yd-1e-7;
colormap(parula);xlim([Interp_xmin,Interp_xmax]); ylim([Interp_ymin Interp_ymax]);
X_meshgrid = reshape(node1(:,1),nx,ny);
Y_meshgrid = reshape(node1(:,2),nx,ny);
U_meshgrid = reshape(p,nx,ny);
hold on;
contour(X_meshgrid,Y_meshgrid,U_meshgrid')
view([50 30]);
filename = ['Curve_',num2str(nnn),'_solution_1'];
saveas(gcf, filename,'epsc');


%slice
[x_rect, y_rect] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
X_meshgrid = x_rect';
Y_meshgrid = y_rect';
U_meshgrid = reshape(p,nx,ny);

filename1 = ['Curve_',num2str(nnn),'_1'];
save(filename1);
[X_slice1, Y_slice1] = meshgrid(linspace(xa,xb,1000),0.5);
R_slice1 = sqrt((X_slice1 - X_slice1(1)).^2+(Y_slice1-Y_slice1(1)).^2);
U_slice1 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice1,Y_slice1,'linear');
save(filename1,'X_slice1','Y_slice1','R_slice1','U_slice1','-append') 
figure;plot(R_slice1,U_slice1,'LineWidth',1.5);

[X_slice2, Y_slice2] = meshgrid(0.4, linspace(yc,yd,1000));
R_slice2 = sqrt((X_slice2 - X_slice2(1)).^2+(Y_slice2-Y_slice2(1)).^2);
U_slice2 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice2,Y_slice2,'linear');
save(filename1,'X_slice2','Y_slice2','R_slice2','U_slice2','-append') 
figure;plot(R_slice2,U_slice2,'LineWidth',1.5);

% t = linspace(0, 1, 1000);
% X_line = t;
% Y_line = t;
% R_line = sqrt(X_line.^2 + Y_line.^2);
% U_line = interp2(X_meshgrid', Y_meshgrid', U_meshgrid', X_line, Y_line, 'linear');

toc


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