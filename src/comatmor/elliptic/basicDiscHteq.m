% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 20.02.2015
% Linked to model heatequation.m


%model = ModelUtil.model('Model2');
% Create function handles
% Usage in command prompt: 
% heatDisc = heatDisc
% solutions = heatDisc.compute(...)
% heatDisc.visualize(...)

function  basicDiscHeateq = basicDiscHeateq()

    basicDiscHeateq.compute = @compute;
    basicDiscHeateq.visualize = @visualize;
end

function solutions= compute(model)

% initialize right python interpreter
pyversion /home/310191226/pymorDir/virt/bin/python

% Get basic modelinfo
modelinfo = mphmodel(model)

% parameters
global sol pg;
% fix solver nod
sol = 'sol1';
% plotgroup
pg = 'pg1';
% set num_samples(the reduced basis generation should consider)
num_samples = 10;

% names of to varying parameters
parameter = 'c';

% user defined parameter_ranges for parameterfile
paramRanges = ['(1,50)';'(1,50)'];

% set parameter_set and write it to disc
global parameter_set

parameter_set = [[1.0,40.0];[40.0,1.0]];
save('parameter_set.mat','parameter_set')
parameterName = '"parameter_set"';
parameterPath = '"parameter_set.mat"';


% Deactivate internal dofs for simple ellip problem
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% AFFINE DECOMPOSITION
modelPhysics = model.physics(modelinfo.physics);
% Get dependent parameters
%mphgetproperties(modelPhysics.feature('hteq1'));

% set variables to 0 and get matrices 
% use try/error block to get all involved equations
% stiffness matrix Kc + load vectors Lc(RHS)
numb = 0;
for i=1:10
    try
        str = ['hteq',int2str(i)];
        modelPhysics.feature(str);
    catch
        fprintf('There are %d heatequations given.\n',(i-1))
        numb = i-1;
        break
    end
end

% save matrixNames and filenames (path) to create parameterfile
matrixNames = '';
matrixPaths = '';
% Save matrices
for i=1:numb
    str = ['hteq',int2str(i)];
    modelPhysics.feature(str).set(parameter,1);
    for j=1:numb
        if j~= i
            strj = ['hteq',int2str(j)];
            modelPhysics.feature(strj).set(parameter,0);
        end
    end
    % Get matrices from comsol
    MA = mphmatrix(model ,sol, 'Out', {'Kc','Lc'},'initmethod','init');
    % Save matrices to harddisc
    KName = ['Kc',int2str(i)];
    LName = ['Lc',int2str(i)];
    matrixNames = [matrixNames; ['"',KName,'"']; ['"',LName,'"']];
    Kc = matlab.lang.makeValidName(KName);
    Lc = matlab.lang.makeValidName(LName);
    eval([Kc '= MA.Kc;']);
    eval([Lc '= MA.Lc;']);
    save([KName,'.mat'],KName)
    save([LName,'.mat'],LName)
    matrixPaths=[matrixPaths; ['"',KName,'.mat"']; ['"',LName,'.mat"']];
    % delete matrices
    eval([Kc '= [];']);
    eval([Lc '= [];']);
    Kc = [];
    Lc = [];
end

% Go to default (later save state before perhaps?)
modelPhysics.feature('hteq1').set(parameter,1);
modelPhysics.feature('hteq2').set(parameter,1);

% Get other components
MA = mphmatrix(model ,sol, ...
'Out', {'Null','ud','uscale'},...
'initmethod','init');

% Create parameterfile for given problem
paramFile = fopen('parameterHeateq.py','w');
fprintf(paramFile,'matfile = {');
% for matDict
length(matrixNames(:,1));

for i=1:length(matrixNames(:,1))
    % Insert correct paramRanges (parameter are the same for 
    j = round(i/2);
    % create ci for parameternames in parameterType in pymor
    varName = ['"c',int2str(j),'"'];
    fprintf(paramFile,[matrixNames(i,:),':','(',matrixPaths(i,:),',','[',varName,']',',','[1]',',','[',paramRanges(j,:),'])']);
    if i==length(matrixNames(:,1))
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,'}\n');
% for stiffNames
fprintf(paramFile,'stiffNames=(');
for i=1:2:length(matrixNames(:,1))
    fprintf(paramFile,matrixNames(i,:));
    if i==length(matrixNames(:,1))-1
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,')\n');
% for rhsNames
fprintf(paramFile,'rhsNames=(');
for i=2:2:length(matrixNames(:,1))
    fprintf(paramFile,matrixNames(i,:));
    if i==length(matrixNames(:,1))
        break
    end
    fprintf(paramFile,',');
end
fprintf(paramFile,')\n');
%parameterSetfile
fprintf(paramFile,['parameterSetfile={',parameterName,':',parameterPath,'}']);


% Call python script
samples = ['--samples=',int2str(num_samples)];
system(['source /home/310191226/pymorDir/virt/bin/activate && python startEllipticRBHeateq.py',' ',samples])

% Load solutions from harddisk as struct M, ensure right numbering
for i=1:length(parameter_set(:,1))
    name=['mu',int2str(i)];
    inter = load('RBsolutions.mat',name);
    M.(name) = inter.(name);
end


% Calculate final solution(scale+boundary conditions)
names = fieldnames(M);

% introduce new variable solutions to fit (possibly) adjustet size
for i=1:numel(names)
    M.(names{i}) = M.(names{i})';
    solutions.(names{i})= zeros(length(MA.Null),numel(names));
end

for i=1:numel(names)
    solutions.(names{i})(:,j)=MA.Null*M.(names{i})(:,j);
    solutions.(names{i})(:,j)=solutions.(names{i})(:,j)+MA.ud;
    solutions.(names{i})(:,j)=(1+solutions.(names{i})(:,j)).*MA.uscale;
    % There was (1+...) in last eq 
end

end

% Set and visualize solution in comsol and matlab
function visualize(solutions,model,sel)

% get global variables
global sol pg parameter_set;

names = fieldnames(solutions);

model.sol(sol).setU(solutions.(names{sel}));
model.sol(sol).createSolution;
mphplot(model,pg,'rangenum',1);
% Give figure title of seen parameter
caption = ['Solution for parameter = (',num2str(parameter_set(sel,:)),')'];
title(caption)
end
