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
% fix solver nod
sol = 'sol2';
% plotgroup
pg = 'pg2';
% mphsearch

% Deactivate internal dofs for simple ellip problem
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% Get initial solution
u0 = mphgetu(model,'solnum',1);
save('u0.mat','u0')

% AFFINE DECOMPOSITION
modelPhysics = model.physics(modelinfo.physics);
% Get dependent parameters
%mphgetproperties(modelPhysics.feature('hteq1'));

% set variables to 0 and get matrices 
% nam = ['hteq',int2str(i)];
% stiffness matrix Kc + load vectors Lc(RHS)
modelPhysics.feature('hteq1').set('c',1);
modelPhysics.feature('hteq2').set('c',0);

MA = mphmatrix(model ,sol, ...
'Out', {'Kc','Lc'},'initmethod','init');

Kc1 = MA.Kc;
Lc1 = MA.Lc;

modelPhysics.feature('hteq1').set('c',0);
modelPhysics.feature('hteq2').set('c',1);

MA = mphmatrix(model ,sol, ...
'Out', {'Kc','Lc'},'initmethod','init');

Kc2 = MA.Kc;
Lc2 = MA.Lc;

% Go to default (later save state before perhaps?)
modelPhysics.feature('hteq1') .set('c',1);
modelPhysics.feature('hteq2').set('c',1);

% Get other components
MA = mphmatrix(model ,sol, ...
'Out', {'Dc','Null','ud','uscale'},...
'initmethod','init');

% Damping matrix
Dc = MA.Dc;
% Save matrices to harddisk as .mat file 
save('Kc1.mat','Kc1')
save('Kc2.mat','Kc2')
save('Lc1.mat','Lc1')
save('Lc2.mat','Lc2')
save('Dc.mat','Dc')

% Call python script
system('source /home/310191226/pymorDir/virt/bin/activate && python startHeatRB.py')

% Load solutions from harddisk 
% As struct M
M = load('RBsolutions.mat');

% Calculate final solution(scale+boundary conditions)
names = fieldnames(M);
steps = length(M.(names{1})(:,1));
for i=1:numel(names)
    M.(names{i}) = M.(names{i})';
end
for i=1:numel(names)
    for j=1:steps
        M.(names{i})(:,j)=MA.Null*M.(names{i})(:,j);
        M.(names{i})(:,j)=M.(names{i})(:,j)+MA.ud;
        M.(names{i})(:,j)=(M.(names{i})(:,j)).*MA.uscale;
        % There was (1+...) in last eq 
    end
end

% Set and visualize solution in comsol and matlab
names{1}
t = model.sol(sol).getPVals
model.sol(sol).setPVals(t)
for i=1:length(t)
    model.sol(sol).setU(i,M.(names{1})(:,i));
end

model.sol(sol).createSolution;
mphplot(model,pg,'rangenum',1);