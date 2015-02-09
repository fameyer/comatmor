% COMSOL-MATLAB-PYMOR interface
% Falk Meyer, 06.02.2015

% Automatize setup of COMSOL-MATLAB Server
% ???

% Get model from COMSOL server
model = ModelUtil.model('Model2'); % Model3 just for example

% Get stiffnessmatrix as 1x1 cell struct
st = mphmatrix(model,'sol1','out',{'K'});

% Get stiffnessmatrix as nxn sparse matrix
S = st.K;

% Save matrix to harddisk as .mat file 
save('stiffnessMatrix.mat','S')

% Get nonzeros indices for matrix entries
% PROBABLY unnecessary
[row,col,data] = find(S);

% Call pymor Routine to convert S into a NumpyMatrixOperator
