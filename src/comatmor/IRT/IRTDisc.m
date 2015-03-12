% COMSOL-MATLAB-PYMOR interface for disc-based communication
% Falk Meyer, 01.03.2015
% Linked to model heatequationTime.m

%model = ModelUtil.model('Model2');
% Create function handles
% Usage in command prompt: 
% heatDisc = heatDisc
% solutions = heatDisc.compute(...)
% heatDisc.visualize(...)

function IRTDisc = IRTDisc()

    IRTDisc.compute = @compute;
    IRTDisc.visualize = @visualize;
end

function solutions= compute(model)

% initialize right python interpreter
pyversion /home/310191226/pymorDir/virt/bin/python

% Get basic modelinfo
modelinfo = mphmodel(model)

% parameters
global sol pg T steps;
% fix solver nod
sol = 'sol1';
% plotgroup
pg = 'pg1';
% set num_samples(the reduced basis generation should consider), endtime and stepsize
T = 5;
steps = 20;
num_samples = 4;

% names of to varying parameters
parameterStiff = 'c';
parameterMass = 'da';

% user defined parameter_ranges for parameterfile
kappaSample_Range = '(1,50)';
rhocSample_Range = '(1,50)';
kappaRest_Range = '(1,1)';
rhocRest_Range = '(1,1)';

paramRanges = [rhocSample_Range;kappaSample_Range];

% set parameter_set and write it to disc
global parameter_set

parameter_set = [[40.0,20.0,1.0]]; % Distinction by semicolon!
save('parameter_set.mat','parameter_set')
parameterName = '"parameter_set"';
parameterPath = '"parameter_set.mat"';

% Deactivate internal dofs to enable comparable results
Shape = model.physics(modelinfo.physics).prop('ShapeProperty');
%Shape.set('boundaryFlux_temperature', 1, '0'); % for ht model
Shape.set('boundaryFlux', 1, '0');

% Read index list of dirichlet indices (has tag index)
indexDirichlet = load('dirichletIndex.mat');
index = indexDirichlet.index;

% % Get initial solution
u0 = mphgetu(model,'solnum',1);
%MA = mphmatrix(model ,sol, 'Out', {'Null'},'initmethod','init');
%u0 = MA.Null'*u0;
% Don't know why i didn't come
%u0 = u0-1;

% Ensure Dirichlet conditions
% for j=1:length(index)
%       i = index(j);
% %      u0(i+1:end+1) = u0(i:end);
%       u0(i) = 50.0;
% end

%u0=20.0*ones(10284,1);
save('u0.mat','u0');
% % strings to store names
% u0Name = '"u0"';
% u0Path = '"u0.mat"';

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
KcSample = MA.K;
LcSample = MA.L;
DcSample = MA.D;

% Ensure Dirichlet conditions
for j=1:length(index)
     i = index(j);
%     %KcSample(i+1:end+1,:) = KcSample(i:end,:);
      KcSample(i,:) = 0;
%     %KcSample(:,i+1:end+1) = KcSample(:,i:end);
      %KcSample(:,i) = 0;
%     KcSample(i,i) = 1;
%     %DcSample(i+1:end+1,:) = DcSample(i:end,:);
     DcSample(i,:) = 0;
%     %DcSample(:,i+1:end+1) = DcSample(:,i:end);
%     %DcSample(:,i) = 0;
%     DcSample(i,i) = 1;
%     %LcSample(i+1:end+1) = LcSample(i:end);
     LcSample(i) = 0.0;
 end

save('KcSample.mat', 'KcSample');
save('LcSample.mat', 'LcSample');
save('DcSample.mat', 'DcSample');

% Get matrices from the rest
modelPhysics.feature('hteq2').set(parameterStiff,0);
modelPhysics.feature('hteq2').set(parameterMass,0);

modelPhysics.feature('hteq1').set(parameterStiff,c1);
modelPhysics.feature('hteq1').set(parameterMass,da1);
modelPhysics.feature('hteq3').set(parameterStiff,c3);
modelPhysics.feature('hteq3').set(parameterMass,da3);
modelPhysics.feature('hteq4').set(parameterStiff,c4);
modelPhysics.feature('hteq4').set(parameterMass,da4);

MA = mphmatrix(model ,sol, 'Out', {'K','L','D'},'initmethod','init');
Kc = MA.K;
Lc = MA.L;
Dc = MA.D;

% Ensure Dirichlet conditions
 for j=1:length(index)
     i = index(j);
%     %Kc(i+1:end+1,:) = Kc(i:end,:);
     Kc(i,:) = 0;
%     %Kc(:,i+1:end+1) = Kc(:,i:end);
%     %Kc(:,i) = 0;
     Kc(i,i) = 1;
%     %Dc(i+1:end+1,:) = Dc(i:end,:);
     Dc(i,:) = 0;
%     %Dc(:,i+1:end+1) = Dc(:,i:end);
%     %Dc(:,i) = 0;
%     Dc(i,i) = 0;
%     %Lc(i+1:end+1) = Lc(i:end);
     Lc(i) = 0.0;
 end

save('Kc.mat', 'Kc');
save('Lc.mat', 'Lc');
save('Dc.mat', 'Dc');

% Go to default values
modelPhysics.feature('hteq2').set(parameterStiff,1);
modelPhysics.feature('hteq2').set(parameterMass,1);

% Get other components
MA = mphmatrix(model ,sol, ...
'Out', {'Null','ud','uscale'},...
'initmethod','init');

% Call python script
endtime =['--endtime=',num2str(T)];
step_number = ['--steps=',int2str(steps)];
samples = ['--samples=',int2str(num_samples)];
system(['source /home/310191226/pymorDir/virt/bin/activate && python startIRTRB.py',' ',endtime,' ',step_number,' ',samples]);

% Load solutions from harddisk as struct M, ensure right numbering
% As struct M
for i=1:length(parameter_set(:,1))
    name=['mu',int2str(i)];
    inter = load('RBsolutions.mat',name);
    M.(name) = inter.(name);
end

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
        solutions.(names{i})(:,j)=M.(names{i})(:,j);%MA.Null*M.(names{i})(:,j);
        %solutions.(names{i})(:,j)=solutions.(names{i})(:,j)+MA.ud;
        %solutions.(names{i})(:,j)=(1+solutions.(names{i})(:,j)).*MA.uscale;
        % There was (1+...) in last eq 
        %solutions.(names{i})(:,j)=20+solutions.(names{i})(:,j);
    end
end

end

% Set and visualize solution in comsol and matlab
% sel is the index of the solution you want to visualize
function visualize(model,solutions,sel)

% get global variables
global sol pg steps T parameter_set;

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
caption = ['Solution for parameter = (',num2str(parameter_set(sel,:)),')'];
title(caption)
end