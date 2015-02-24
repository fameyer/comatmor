% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 20.02.2015
% Linked to model heatequation.m

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
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% AFFINE DECOMPOSITION
modelPhysics = model.physics(modelinfo.physics);
% Get dependent parameters
%mphgetproperties(modelPhysics.feature('hteq1'));

% set variables to 0 and get matrices 
% nam = ['hteq',int2str(i)];
% first stiffness matrix
modelPhysics.feature('hteq1').set('c',1);
modelPhysics.feature('hteq2').set('c',0);

MA = mphmatrix(model ,'sol1', ...
'Out', {'Kc'},'initmethod','init');

Kc1 = MA.Kc;

% second stiffness matrix
modelPhysics.feature('hteq1').set('c',0);
modelPhysics.feature('hteq2').set('c',1);

MA = mphmatrix(model ,'sol1', ...
'Out', {'Kc'},'initmethod','init');

Kc2 = MA.Kc;

% Go to default (later save state before perhaps?)
modelPhysics.feature('hteq1') .set('c',1);
modelPhysics.feature('hteq2').set('c',1);

% Get other components
MA = mphmatrix(model ,'sol1', ...
'Out', {'Lc','Null','ud','uscale'},...
'initmethod','init');

% Get right-hand side 
Lc = MA.Lc;

% Save matrices to harddisk as .mat file 
save('Kc1.mat','Kc1')
save('Kc2.mat','Kc2')
save('rhsHeateq.mat','Lc')

% Call python script
system('source /home/310191226/pymorDir/virt/bin/activate && python startEllipticRBHeateq.py')

% Load solutions from harddisk 
% As struct M
M = load('RBsolutions.mat');

% Calculate final solution(scale+boundary conditions)
names = fieldnames(M);
for i=1:numel(names)
    M.(names{i})=MA.Null*M.(names{i})';
    M.(names{i})=M.(names{i})+MA.ud;
    M.(names{i})=(M.(names{i})).*MA.uscale;
    % There was (1+...) in last eq ?
end

% Set and visualize solution in comsol and matlab
names{1}
model.sol('sol1').setU(M.(names{1}));
model.sol('sol1').createSolution;
mphplot(model,'pg1','rangenum',1);