% set the seed to 0
rng(0);
clear;

% Load model
mdlLoader = iDynTree.ModelLoader();
params.modelFile = 'twoLinks.urdf';
mdlLoader.loadModelFromFile(params.modelFile);

% Create berdyHelper 
% To simulate the classical case, disable external forces and use joint 
% torques and accelerometers as sensors 
berdyOptions = iDynTree.BerdyOptions();
berdyOptions.includeAllNetExternalWrenchesAsSensors          = false;
berdyOptions.includeAllNetExternalWrenchesAsDynamicVariables = false;
berdyOptions.includeAllJointAccelerationsAsSensors           = true;
berdyOptions.includeAllJointTorquesAsSensors                 = true;
berdy = iDynTree.BerdyHelper();
berdy.init(mdlLoader.model(),mdlLoader.sensors(),berdyOptions);

% Create dynComp 
dynComp = iDynTree.DynamicsComputations();
dynComp.loadRobotModelFromFile(params.modelFile);

% Create Buffers of iDynTree objects   
buffers.jntPos = iDynTree.JointPosDoubleArray(berdy.model());
buffers.jntVel = iDynTree.JointDOFsDoubleArray(berdy.model());
buffers.jntAcc = iDynTree.JointDOFsDoubleArray(berdy.model());
buffers.jntTrqs = iDynTree.JointDOFsDoubleArray(berdy.model());
buffers.baseReactionForce = iDynTree.Wrench();
buffers.jntPos.zero();
buffers.jntVel.zero();
buffers.jntAcc.zero();
buffers.jntTrqs.zero();
buffers.fullDynRegrId = iDynTree.MatrixDynSize();
buffers.fullCadParams = iDynTree.VectorDynSize();
buffers.fullEstimatedParams = iDynTree.VectorDynSize();

% Check that dynComp (old KDL-based code) and model serialization of the
% inertial parameters are coincident. In general this is not true (because 
% the link serialization depends on the parser), but for simple chain
% models used in this script this should hold
% 
dynComp.getModelDynamicsParameters(buffers.fullCadParams);
CADParamsDynComp = buffers.fullCadParams.toMatlab();
berdy.model().getInertialParameters(buffers.fullCadParams);
CADParamsBerdy = buffers.fullCadParams.toMatlab();
assert(norm(CADParamsDynComp-CADParamsBerdy) < 1e-10);


% Gravity 
gravM = [0.0 0.0 -10.0];
buffers.grav = iDynTree.Vector3();
buffers.grav.fromMatlab(gravM);
buffers.gravSpatialAcc = iDynTree.SpatialAcc();
buffers.gravSpatialAcc.fromMatlab([gravM zeros(1,3)]);


% Update kinematics 
berdy.updateKinematicsFromTraversalFixedBase(buffers.jntPos,buffers.jntVel,buffers.grav);

% Get matrices 
buffers.Yid = iDynTree.MatrixDynSize();
buffers.Did = iDynTree.MatrixDynSize();
buffers.bYid = iDynTree.VectorDynSize();
buffers.bDid = iDynTree.VectorDynSize();
berdy.resizeAndZeroBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
berdy.getBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
D  = buffers.Did.toMatlab();
bD = buffers.bDid.toMatlab();
Y  = buffers.Yid.toMatlab();
bY = buffers.bYid.toMatlab();

nrOfBerdyDynEquations = size(bD,1);
[nrOfBerdySensors,nrOfBerdyDynVariables] = size(Y);

% Generate training trajectory 
% parameters 
dofs = berdy.model().getNrOfDOFs();
params.freq = 2.0; 
params.A    = 1.0;
params.w    = repmat(2*pi*params.freq,dofs,1);
params.ddqRelNoise = 0.1;
params.trqsRelNoise = 0.1;
params.stdDevDynEq  = 1e-3;
params.stdDevDynVariablesPrior = 1e3;

dofs = berdy.model().getNrOfDOFs();
dataset.t = 0:0.02:3;

dataset.q = sin(params.w*dataset.t);
dataset.dq = params.w*cos(params.w*dataset.t);
dataset.ddq = -params.w*params.w*sin(params.w*dataset.t);

[dofs,nrOfSamples] = size(dataset.q);

% Compute dataset.trq throught RNEA
dataset.trqs = zeros(size(dataset.dq));
for i = 1:nrOfSamples
    buffers.jntPos.fromMatlab(dataset.q(:,i));
    buffers.jntVel.fromMatlab(dataset.dq(:,i));
    buffers.jntAcc.fromMatlab(dataset.ddq(:,i));
    
    dynComp.setRobotState(buffers.jntPos,buffers.jntVel,buffers.jntAcc,buffers.gravSpatialAcc);
    dynComp.inverseDynamics(buffers.jntTrqs,buffers.baseReactionForce);
    
    dataset.trqs(:,i) = buffers.jntTrqs.toMatlab();
end

% Generate simulated measurements, for now just using jnt torques and joint
% accelerations
measurements.ddq = zeros(size(dataset.ddq));
measurements.trqs = zeros(size(dataset.trqs));

params.stdDev_ddq = params.ddqRelNoise*std(dataset.ddq')';
params.stdDev_trqs = params.trqsRelNoise*std(dataset.trqs')';
for i = 1:nrOfSamples
    measurements.ddq(:,i) = dataset.ddq(:,i) + normrnd(0,params.stdDev_ddq,dofs,1);
    measurements.trqs(:,i) = dataset.trqs(:,i) + normrnd(0,params.stdDev_trqs,dofs,1);
end

paramsStr = sprintf('ddq %.2e trqs %.2e dynEq %.2e dynVar %.2e',params.stdDev_ddq,params.stdDev_trqs,params.stdDevDynEq,params.stdDevDynVariablesPrior);

%% Compute identifiable parameters from the specified dataset (currently hardcode to always give just the more relevant parameter)
identifibleParamsMatrix = computeJointTorquesExperimentalIdentifiableParams(dynComp,buffers,params,dataset,measurements);
[nrOfIdentifiableParams,nrOfTotalParams] = size(identifibleParamsMatrix);
% Currently we are assuming to just identify just one inertial parameter 
assert(nrOfIdentifiableParams == 1);
dynComp.getModelDynamicsParameters(buffers.fullCadParams);
groundTruthValues = identifibleParamsMatrix*buffers.fullCadParams.toMatlab();

% For simplifyng the log data exploration, we make sure that the selected
% parameters is positive 
if( groundTruthValues < 0.0 )
    identifibleParamsMatrix = -identifibleParamsMatrix;
end
groundTruthValues = identifibleParamsMatrix*buffers.fullCadParams.toMatlab();


%% Estimate inertial parameters, return cell array of estimated parameters for number of used samples 
classicalEstimationResults = getClassicalEstimationResults(dynComp,buffers,dataset,measurements,identifibleParamsMatrix);

%% Plot error in classical estimation results 
classicalEstimationOnRelevantParameters = classicalEstimationResults;

classicalEstimationRelativeErrors = diag(1./groundTruthValues)*(classicalEstimationOnRelevantParameters-repmat(groundTruthValues,1,nrOfSamples));
% Drop the first few samples 
figure;
title('Relative error on inertial parameters identification')
plot(dataset.t(3:end),classicalEstimationRelativeErrors(:,3:end)');
hold on;
plot(dataset.t(3:end),zeros(size(dataset.t(3:end))));
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',paramsStr,'FitBoxToText','on');
title('Relative error on inertial parameters identification')
% ylim([-3,3]);

%% Load covariances matrices 
covs.cov_e_given_d = params.stdDevDynEq*params.stdDevDynEq*eye(nrOfBerdyDynEquations,nrOfBerdyDynEquations);
covs.cov_d         = params.stdDevDynVariablesPrior*params.stdDevDynVariablesPrior*eye(nrOfBerdyDynVariables,nrOfBerdyDynVariables);
% Sensors should be 2*dofs : first the accelerations then the torques 
% This should be modified to have a nice introspection mechanism, but for
% now this is ok 
assert(nrOfBerdySensors == 2*berdy.model().getNrOfDOFs());
nrOfDOFs = berdy.model().getNrOfDOFs();
cov_ddq = params.stdDev_ddq*params.stdDev_ddq*eye(nrOfDOFs,nrOfDOFs);
cov_trqs = params.stdDev_ddq*params.stdDev_trqs*eye(nrOfDOFs,nrOfDOFs);
covs.cov_y_given_d   = blkdiag(cov_ddq,cov_trqs);
measurements.y = [measurements.ddq;measurements.trqs];

%% Let's define the loglikeliwood
ll = @(pi) getLogLikeLihoodOfInertialParameters(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,pi);
minusll = @(pi) -ll(pi);

%% Plot log likeliwood of estimation 
llClassicalEstimation = zeros(nrOfSamples,1);
lClassicalEstimation  = zeros(nrOfSamples,1);
for sampleIdx = 1:nrOfSamples
    fprintf('Computing ll for sample cad estimate at sample %d of of %d\n',sampleIdx,nrOfSamples)
    [llClassicalEstimation(sampleIdx),lClassicalEstimation(sampleIdx)] = ll(classicalEstimationResults(:,sampleIdx));
end

% Get ll of totally wrong data 
llGroundTruth = getLogLikeLihoodOfInertialParameters(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,groundTruthValues);

figure;
plot(dataset.t(3:end),llClassicalEstimation(3:end));
hold on;
plot(dataset.t(3:end),llGroundTruth*ones(size(dataset.t(3:end))));
annotation('textbox',dim,'String',paramsStr,'FitBoxToText','on');
title('Evolution of LogLikelihood of inertial parameters estimated');

% Get ll of actual data 
% llGroundTruth = getLogLikeLihoodOfInertialParameters(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,groundTruthValues);

% Let's explore the face of the ll function 
inPar = 1.07:0.01:2;
llExploration = zeros(size(inPar));
for sampleIdx = 1:size(inPar,2);
    fprintf('Computing ll for inertial parameters %f at sample %d of of %d\n',inPar(sampleIdx),sampleIdx,size(inPar,2))
    llExploration(sampleIdx) = ll(inPar(sampleIdx));
end

% Let's find the maximum of the ll function
inertialParamEstimatedMinizingLL = fminsearch(minusll,groundTruthValues);

figure;
plot(inPar,(llExploration),'b.');
hold on;
plot([groundTruthValues groundTruthValues], [min((llExploration)) max((llExploration))], 'g');
plot([classicalEstimationResults(:,nrOfSamples) classicalEstimationResults(:,nrOfSamples)], [min((llExploration)) max((llExploration))], 'r');
plot([inertialParamEstimatedMinizingLL inertialParamEstimatedMinizingLL], [min((llExploration)) max((llExploration))], 'b');
annotation('textbox',dim,'String',paramsStr,'FitBoxToText','on');
title('LogLikelihood function');

