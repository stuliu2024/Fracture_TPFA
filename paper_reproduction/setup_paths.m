function setup_paths
root=fileparts(mfilename('fullpath'));
mrst=fullfile(root,'dependencies','mrst-2023a');
parts={'gridprocessing','params','plotting','solvers','utils'};
for k=1:numel(parts)
    folders=strsplit(genpath(fullfile(mrst,parts{k})),pathsep);
    for j=1:numel(folders)
        if ~isempty(folders{j}) && isempty(strfind(folders{j},'3rdparty')) && isempty(strfind(folders{j},'octave_only'))
            addpath(folders{j});
        end
    end
end
mrstPath('reregister','incomp',fullfile(mrst,'modules','incomp'), ...
    'linearsolvers',fullfile(mrst,'modules','linearsolvers'));
mrstModule add incomp linearsolvers;
addpath(fullfile(mrst,'utils','3rdparty','agmg','Matlab'));
addpath(fullfile(root,'dependencies','quadrature'));
gravity reset off;
end
