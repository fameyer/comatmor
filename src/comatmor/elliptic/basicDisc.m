% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 06.02.2015

% Automatize setup of COMSOL-MATLAB Server - write function
% ???

% initialize right python interpreter
pyversion /home/310191226/pymorDir/virt/bin/python

% Get model from COMSOL server
model = ModelUtil.model('Model2'); % Model2 just for example
% Get basic modelinfo
modelinfo = mphmodel(model)

% mphsearch

% Deactivate internal dofs for simple ellip problem
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model

MA = mphmatrix(model ,'sol1', ...
'Out', {'Kc','Lc','Null','ud','uscale'},...
'initmethod','init');

% Get stiffnessmatrix as nxn sparse matrix
S = MA.Kc;

% Get right-hand side 
L = MA.Lc;

% Save matrices to harddisk as .mat file 
save('matrix.mat','S')
save('rhs.mat','L')

% Call python script
system('source /home/310191226/pymorDir/virt/bin/activate && python startEllipticRB.py')

% Load solutions from harddisk 
% As struct M
M = load('RBsolutions.mat');

% Calculate final solution(scale)
% Shortens to statement below due to MA.Null = 1 and MA.ud = 0 for THIS
% case
% Uc = MA.Kc\MA.Lc;
%U1 = (1+Uc).*MA.uscale;
names = fieldnames(M);
for i=1:numel(names)
    M.(names{i})=(1+M.(names{i})').*MA.uscale;
end

% Set and visualize solution in comsol and matlab
model.sol('sol1').setU(M.(names{1}));
model.sol('sol1').createSolution;
mphplot(model,'pg1','rangenum',1);