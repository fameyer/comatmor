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

% Get parameters and variables
%modelparaInfo = mphgetexpressions(model)
% get/change physics parameter
% ht=model.physics('ht')
% mphgetproperties(ht.feature('solid1'))
% model.param.set('k',2000)
% mphgetexpressions(model.param)
% print mesh
% mphmesh(model)

% Get stiffnessmatrix as 1x1 cell struct
st = mphmatrix(model,'sol1','out',{'K'});

% Get stiffnessmatrix as nxn sparse matrix
S = st.K;

% Get right-hand side 

% Get solution of model
%U = mphgetu(model);

% Get RHS?
%lt = mphmatrix(model,'sol1','out',{'L'});

% Call matlab pymor Interface
% Save matrix to harddisk as .mat file 
save('stiffnessMatrix.mat','S')

% Find possibility to set parameters

% Call python script
system('source /home/310191226/pymorDir/virt/bin/activate && python startRB.py')
% Change conductivity of model prob
%model.material('mat1').propertyGroup('def').set('thermalconductivity', '2000');

% Load solutions from harddisk
load('RBsolutions.mat')

% Add solution to comsol model (does this work?)
%model.sol().create(NAME)