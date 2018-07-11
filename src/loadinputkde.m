function inputkde = loadinputkde() %#ok<STOUT>
%%% loads the input kde object
%% get path
% first to the data_import script then back out of `src` 
% to the data folder
pathstr = fileparts(which('data_import'));
pathstr = fullfile(pathstr(1:end-4),'data');

%% load input kde object and sample
load(fullfile(pathstr,'inputkde.mat'),'inputkde');