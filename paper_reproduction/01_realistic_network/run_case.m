clear; close all;
root=fileparts(mfilename('fullpath'));addpath(root);
nx=175; ny=150;
solve_realcase(nx,ny);
plot_realcase;
make_efficiency_table;
