clear;close all;
root=fileparts(mfilename('fullpath'));addpath(root);
% Main text: 20,35,50. Supporting Information: 9,17,35.
for N=[9 17 20 35 50]
    run_regular_tpfa(N);
end
plot_regular_network;
plot_si_profiles;
