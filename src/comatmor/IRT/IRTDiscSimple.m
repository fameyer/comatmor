% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 01.03.2015
% Linked to model IRTsimple.mph

% model = ModelUtil.model('Model2');

% Usage in command prompt: 
% IRTDisc = IRTDiscSimple
% solutions = IRTDisc.compute(...)
% IRTDisc.visualize(...)

function IRTDisc = IRTDiscSimple()

    IRTDisc.compute = @compute;
    IRTDisc.visualize = @visualize;
end

function solutions= compute(model, parameter_set, num_samples, max_extensions, target_error)

% Get basic modelinfo
modelinfo = mphmodel(model)

% parameters
global sol T steps;
% fix solver node
sol = 'sol1';
% names of to varying parameters
parameterStiff = 'c';
parameterMass = 'da';

% user defined parameter_ranges for parameterfile
daSample_Range = '(1,50)';
cSample_Range = '(1,50)';

% set num_samples(the reduced basis generation should consider), endtime T
% and number of steps
T = 5.0;
steps = 20;

if (~exist('num_samples','var'))
    num_samples = 4;
end
if (~exist('max_extensions','var'))
    max_extensions = 10;
end
if (~exist('target_error','var'))
    target_error = 1e-10;
end

% set parameter_set and write it to disc
% parameter structure: [c,da,k] k==1.0 
if (~exist('parameter_set','var'))
    parameter_set = [[40.0,20.0,1.0];[1.0,40.0,1.0]]; % Distinction by semicolon!
end

save('parameter_set.mat','parameter_set')
parameterName = '"parameter_set"';
parameterPath = '"parameter_set.mat"';

% Deactivate internal dofs to enable comparable results
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% Get initial solution
u0 = mphgetu(model,'solnum',1);

% build up signature and test if new matrices are needed for computation
sign=['IRT_D_K_L_DSample_LSample_KSample_',int2str(num_samples),'_',int2str(steps),'_',num2str(T,'%1.1f'),'_',int2str(length(u0)),'_',int2str(max_extensions)];
% corresponding indicator
doSave = 1;

% Check if sign was already defined
try
signId = fopen('sign.txt');
tline = fgetl(signId);
    while ischar(tline)
        %disp(tline)
        if strcmp(tline,sign)
            doSave = 0;
            break
        end
        tline = fgetl(signId);
    end
fclose(signId);
catch
    doSave = 1;
end

if doSave == 1
    
disp('Building affine decomposition and saving matrices to harddisc...')

% Save whole system matrix for error_norm
MA = mphmatrix(model ,sol, 'Out', {'K'},'initmethod','init');
K = MA.K;
save('KNorm.mat','K');


% strings to store names
u0Name = '"u0"';
u0Path = '"u0.mat"';

% save initial data
save('u0.mat','u0');

% save dirichlet indices
MA = mphmatrix(model ,sol, 'Out', {'ud'},'initmethod','init');
ud = MA.ud;
index = [];
for i=1:length(ud)
    if ud(i)~=0
        try
            index(end+1)=i;
        catch
            index(1)=i;
        end
     end
 end
save('dirichletIndex.mat','index');

% Read index list of dirichlet indices (has tag index)
%indexDirichlet = load('dirichletIndex.mat');
%index = indexDirichlet.index;

% AFFINE DECOMPOSITION
modelPhysics = model.physics(modelinfo.physics);
% Save previous parameters of non-sample area
m = mphgetproperties(modelPhysics.feature('hteq1'));
c1 = m.c;
da1 = m.da;
m = mphgetproperties(modelPhysics.feature('hteq3'));
c3 = m.c;
da3 = m.da;
m = mphgetproperties(modelPhysics.feature('hteq4'));
c4 = m.c;
da4 = m.da;

% Get matrices from sample
% Sample heateq: hteq2
modelPhysics.feature('hteq2').set(parameterStiff,1);
modelPhysics.feature('hteq2').set(parameterMass,1);

modelPhysics.feature('hteq1').set(parameterStiff,0);
modelPhysics.feature('hteq1').set(parameterMass,0);
modelPhysics.feature('hteq3').set(parameterStiff,0);
modelPhysics.feature('hteq3').set(parameterMass,0);
modelPhysics.feature('hteq4').set(parameterStiff,0);
modelPhysics.feature('hteq4').set(parameterMass,0);

MA = mphmatrix(model ,sol, 'Out', {'K','L','D'},'initmethod','init');
KSample = MA.K;
LSample = MA.L;
DSample = MA.D;

% Extract constant value matrices (which don't vary) for K and L
modelPhysics.feature('hteq2').set(parameterStiff,10);
MA = mphmatrix(model ,sol, 'Out', {'K','L'},'initmethod','init');
KConst = MA.K;
LConst = MA.L;

% Adjust matrices
[iSample, jSample, sSample] = find(KSample);
sConst = nonzeros(KConst);

for i=1:length(sSample)
    if round(sSample(i),12) == round(sConst(i),12)
        KSample(iSample(i),jSample(i)) = 0.0;
    end
end

% Adjust matrices
[iSample, jSample, sSample] = find(LSample);
sConst = nonzeros(LConst);

for i=1:length(sSample)
    if sSample(i) == sConst(i)
        LSample(iSample(i),jSample(i)) = 0.0;
    end
end

% The same for mass matrix
modelPhysics.feature('hteq2').set(parameterMass,10);
MA = mphmatrix(model ,sol, 'Out', {'D'},'initmethod','init');
DConst = MA.D;

% Adjust matrices
[iSample, jSample, sSample] = find(DSample);
sConst = nonzeros(DConst);

for i=1:length(sSample)
    if round(sSample(i),12) == round(sConst(i),12)
        DSample(iSample(i),jSample(i)) = 0.0;
    end
end

% Ensure Dirichlet conditions
for j=1:length(index)
     i = index(j);
     KSample(i,:) = 0;
     DSample(i,:) = 0;
     LSample(i) = 0.0;
end

save('KSample.mat', 'KSample');
save('LSample.mat', 'LSample');
save('DSample.mat', 'DSample');

% Get matrices from the rest
modelPhysics.feature('hteq2').set(parameterStiff,0);
modelPhysics.feature('hteq2').set(parameterMass,0);

modelPhysics.feature('hteq1').set(parameterStiff,c1);
modelPhysics.feature('hteq1').set(parameterMass,da1);
modelPhysics.feature('hteq3').set(parameterStiff,c3);
modelPhysics.feature('hteq3').set(parameterMass,da3);
modelPhysics.feature('hteq4').set(parameterStiff,c4);
modelPhysics.feature('hteq4').set(parameterMass,da4);

MA = mphmatrix(model, sol, 'Out', {'K','L','D'},'initmethod','init');
K = MA.K;
L = MA.L;
D = MA.D;

% Ensure Dirichlet conditions
for j=1:length(index)
     i = index(j);
     K(i,:) = 0;
     K(i,i) = 1;
     D(i,:) = 0;
     L(i) = 0.0;
end

save('K.mat', 'K');
save('L.mat', 'L');
save('D.mat', 'D');

% Save names of matrices and corresponding parameters
matrixNames = {'K','L', 'D', 'KSample', 'LSample', 'DSample'};
paramNames = {'k','k','k','c','c','da'};
matrixPaths = {'K.mat', 'L.mat', 'D.mat', 'KSample.mat', 'LSample.mat', 'DSample.mat'};
paramRanges = {'(1,1)','(1,1)','(1,1)',cSample_Range, cSample_Range, daSample_Range};
stiffNames = {'K', 'KSample'};
massNames =  {'D', 'DSample'};
rhsNames = {'L', 'LSample'};

% Create parameterfile for given problem
paramFile = fopen('parameterIRT.py','w');
fprintf(paramFile,'matfile = {');

% Create correct matDict in parameterfile
for i=1:length(matrixNames)
    fprintf(paramFile,['"',matrixNames{i},'"',':','("',matrixPaths{i},'",','["',paramNames{i},'"],','[1]',',','[',paramRanges{i},'])']);
    if i==length(matrixNames)
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,'}\n');
% for stiffNames
fprintf(paramFile,'stiffNames=('); 
for i=1:length(stiffNames)
    fprintf(paramFile,['"',stiffNames{i},'"']);
    if i==(length(stiffNames))
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,')\n');
% for rhsNames
fprintf(paramFile,'rhsNames=('); 
for i=1:length(rhsNames)
    fprintf(paramFile,['"',rhsNames{i},'"']);
    if i==(length(rhsNames))
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,')\n');
% for massNames
fprintf(paramFile,'massNames=('); 
for i=1:length(massNames)
    fprintf(paramFile,['"',massNames{i},'"']);
    if i==length(massNames)
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,')\n');
%u0file
fprintf(paramFile,['u0file={',u0Name,':',u0Path,'}\n']);
%parameterSetfile
fprintf(paramFile,['parameterSetfile={',parameterName,':',parameterPath,'}']);

% End of matrix writing block
else
    disp('Do not save matrices again...')
end

disp('Calling python script...')

% Call python script
endtime =['--endtime=',num2str(T)];
step_number = ['--steps=',int2str(steps)];
samples = ['--samples=',int2str(num_samples)];
extensions = ['--max_extensions=',int2str(max_extensions)];
error = ['--target_error=',num2str(target_error)];
system(['source /home/310191226/pymorDir/virt/bin/activate && python startIRTRB.py',' ',endtime,' ',step_number,' ',samples,' ',extensions,' ',error]);

disp('Load and convert solutions...')

% Load solutions from harddisk as struct M, ensure right numbering
% As struct M
for i=1:length(parameter_set(:,1))
    name=['mu',int2str(i)];
    inter = load('RBsolutions.mat',name);
    M.(name) = inter.(name)';
end

solutions = M;

end

% Set and visualize solution in comsol and matlab
% sel is the index of the solution you want to visualize
function visualize(model,solutions,sel,time,pg)

% get global variables
global sol steps T;
% check whether plotgroup was set
if (~exist('pg','var'))
        pg = 'pg1';
end
if (~exist('time','var'))
        time = 1;
end

names = fieldnames(solutions);

% set time frame
t= 0:T/steps:T;
%t = model.sol(sol).getPVals
model.sol(sol).setPVals(t)

for i=1:(length(t))
    model.sol(sol).setU(i,solutions.(names{sel})(:,i));
end

disp(['Visualize for plotgroup ',pg,' and timestep ',int2str(time)])

model.sol(sol).createSolution;
mphplot(model,pg,'rangenum',1);

parameter_set = load('parameter_set.mat');

% Give figure title of seen parameter
caption = ['Solution for parameter = (',num2str(parameter_set.parameter_set(sel,:)),')'];
title(caption)
end