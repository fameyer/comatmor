% COMSOL-MATLAB-PYMOR interface
% Falk Meyer, 06.02.2015

% Automatize setup of COMSOL-MATLAB Server
% ???

% initialize right python interpreter
pyversion /home/310191226/pymorDir/virt/bin/python

% Get model from COMSOL server
model = ModelUtil.model('Model2'); % Model3 just for example

% Get stiffnessmatrix as 1x1 cell struct
st = mphmatrix(model,'sol1','out',{'K'});

% Get stiffnessmatrix as nxn sparse matrix
S = st.K;

% Get solution of model
U = mphgetu(model);

% Get loadvector (RHS?)
%lt = mphmatrix(model,'sol1','out',{'L'});
%L = lt.L;

% Call matlab pymor Interface
% Save matrix to harddisk as .mat file 
save('stiffnessMatrix.mat','S')

% Get nonzeros indices for matrix entries
%[row,col,data] = find(S);
% Correct for matlab-python count
%row = row - 1;
%col = col - 1;
% Define stationary RB model
%RB = py.comatmor.stationRB();
% RB.addMatrix('Stiffness matrix', py.numpy.array([row']), py.numpy.array([col']), py.numpy.array([data']), 'Diffusion', 1, py.numpy.array([1,3]), py.numpy.array([1,3]))

% Add solution to comsol model (does this work?)
model.sol().create(NAME)