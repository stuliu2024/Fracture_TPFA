function result=load_regular_result(N)
root=fileparts(mfilename('fullpath'));
file=fullfile(root,'results',sprintf('N%d',N),'rdfm_tpfa_solution.mat');
if ~exist(file,'file'),file=fullfile(root,'saved_results',sprintf('N%d',N),'rdfm_tpfa_solution.mat');end
S=load(file,'result');result=S.result;
end
