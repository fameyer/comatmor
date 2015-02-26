% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 20.02.2015
% Linked to model heatequationTime.m

%model = ModelUtil.model('Model2');
% Create function handles
% Usage in command prompt: 
% heatDisc = heatDisc
% solutions = heatDisc.compute(...)
% heatDisc.visualize(...)

function heatDisc = heatDisc()

    heatDisc.compute = @compute;
    heatDisc.visualize = @visualize;
end

function solutions= compute(model)

% initialize right python interpreter
pyversion /home/310191226/pymorDir/virt/bin/python

% Get model from COMSOL server
%model = ModelUtil.model('Model2'); % Model2 just for example

% Get basic modelinfo
modelinfo = mphmodel(model)

% parameters
global sol pg T steps;
% fix solver nod
sol = 'sol1';
% plotgroup
pg = 'pg1';
% set endtime and stepsize
T = 1;
steps = 10;

% names of to varying parameters
parameter = 'c';

% set training_set and write it to disc
global training_set

training_set = [[1.0,40.0];[40.0,1.0]];
save('training_set.mat','training_set')
% mphsearch

% Deactivate internal dofs to enable comparable results
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% Get initial solution
u0 = mphgetu(model,'solnum',1);
MA = mphmatrix(model ,sol, 'Out', {'Null'},'initmethod','init');
u0 = MA.Null'*u0;
% Don't know why i didn't come
%u0 = u0-1;
save('u0.mat','u0');

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
    Kc = matlab.lang.makeValidName(KName);
    Lc = matlab.lang.makeValidName(LName);
    eval([Kc '= MA.Kc;']);
    eval([Lc '= MA.Lc;']);
    save([KName,'.mat'],KName)
    save([LName,'.mat'],LName)
    % delete matrices
    eval([Kc '= [];']);
    eval([Lc '= [];']);
    Kc = [];
    Lc = [];
end

% Go to default (later save state before perhaps?)
modelPhysics.feature('hteq1') .set('c',1);
modelPhysics.feature('hteq2').set('c',1);

% Get other components
MA = mphmatrix(model ,sol, ...
'Out', {'Dc','Null','ud','uscale'},...
'initmethod','init');

% Damping matrix
Dc = MA.Dc;
save('Dc.mat','Dc')

% Call python script
endtime =['--endtime=',num2str(T)];
step_number = ['--steps=',int2str(steps)];
system(['source /home/310191226/pymorDir/virt/bin/activate && python startHeatRB.py',' ',endtime,' ',step_number])

% Load solutions from harddisk 
% As struct M
M = load('RBsolutions.mat');

% Calculate final solution(scale+boundary conditions)
names = fieldnames(M);
num = length(M.(names{1})(:,1));

% introduce new variable solutions to fit (possibly) adjustet size
for i=1:numel(names)
    M.(names{i}) = M.(names{i})';
    solutions.(names{i})= zeros(length(MA.Null),numel(names));
end

for i=1:numel(names)
    for j=1:num
        solutions.(names{i})(:,j)=MA.Null*M.(names{i})(:,j);
        solutions.(names{i})(:,j)=solutions.(names{i})(:,j)+MA.ud;
        solutions.(names{i})(:,j)=(1+solutions.(names{i})(:,j)).*MA.uscale;
        % There was (1+...) in last eq 
    end
end

end


% Set and visualize solution in comsol and matlab
function visualize(solutions,model,sel)

% get global variables
global sol pg steps T training_set;

names = fieldnames(solutions);

% set time frame
t= 0:T/steps:T;
%t = model.sol(sol).getPVals
model.sol(sol).setPVals(t)

for i=1:(length(t))
    model.sol(sol).setU(i,solutions.(names{sel})(:,i));
end

model.sol(sol).createSolution;

mphplot(model,pg,'rangenum',1);
% Give figure title of seen parameter
caption = ['Solution for parameter = (',num2str(training_set(sel,:)),')'];
title(caption)
end