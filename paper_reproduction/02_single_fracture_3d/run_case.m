clear;close all;
root=fileparts(mfilename('fullpath'));addpath(root);
N=100;
solve_tpfa(N);
solve_cg(N);
plot_single_fracture;
