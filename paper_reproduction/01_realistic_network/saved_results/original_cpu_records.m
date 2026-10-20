%% Data: [Cells, CPU(s), MatrixErr(m), FractureErr(f)]

% --- TPFA
TPFA_cells = [70*16, 70*8, 70*4, 70*2, 70*1];
TPFA_t     = [ 24.876059, 7.377786800000000, 2.917892600000000, 1.475651100000000, 0.815240700000000];

% --- CG
CG_cells = [70*16, 70*8, 70*4, 70*2, 70*1];
CG_t     = [ 43.536576500000002, 11.699574300000000, 3.987735200000000, 1.647722700000000, 0.867741200000000];

% --- LDG
LDG_cells = [70*16, 70*8, 70*4, 70*2, 70*1];
LDG_t     = [ 4.228074398000000e+02, 1.482262497000000e+02, 23.091063299999998, 5.795090600000000, 1.726051500000000];

%% Sort by cells ascending (coarse -> fine)
[TPFA_cells, idx] = sort(TPFA_cells, 'ascend');
TPFA_t  = TPFA_t(idx);  

[CG_cells, idx] = sort(CG_cells, 'ascend');
CG_t  = CG_t(idx);  

[LDG_cells, idx] = sort(LDG_cells, 'ascend');
LDG_t  = LDG_t(idx);  

%% Define mesh size h
% 假设 cells = 每方向网格数 N，且区域长度 L=1，则 h = 1/N
h = 700 ./ TPFA_cells;   % 三组 cells 相同，统一用 TPFA_cells 即可

%% (3) CPU time vs h
figure;
loglog(h, TPFA_t, '-o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
loglog(h, CG_t,   '-s', 'LineWidth', 1.5, 'MarkerSize', 4);
loglog(h, LDG_t,  '-^', 'LineWidth', 1.5, 'MarkerSize', 4);
grid off; box on;
set(gca,'XDir','reverse');
xlabel('h');
ylabel('CPU time (s)');
% title('CPU time vs h');
legend('RDFM-TPFA','RDFM-CG','RDFM-LDG','Location','northwest');
