% set the seed to 0
clear;
rng(2);

% Load model
mdlLoader = iDynTree.ModelLoader();
params.modelFile = 'fourLinks.urdf';
mdlLoader.loadModelFromFile(params.modelFile);

% Create berdyHelper (with just joint torques + joint acceleration as
% sensors)
% To simulate the classical case, disable external forces and use joint 
% torques and accelerometers as sensors 
berdyOptions = iDynTree.BerdyOptions();
berdyOptions.includeAllNetExternalWrenchesAsSensors          = false;
berdyOptions.includeAllNetExternalWrenchesAsDynamicVariables = false;
berdyOptions.includeAllJointAccelerationsAsSensors           = true;
berdyOptions.includeAllJointTorquesAsSensors                 = true;
berdy = iDynTree.BerdyHelper();
berdy.init(mdlLoader.model(),iDynTree.SensorsList(),berdyOptions);

% Create a berdyHelper with joint torques, joint acceleration, 
% link accelerations and one six axis F/T sensor in the base
berdyOptionsRedundantSensors = iDynTree.BerdyOptions();
berdyOptionsRedundantSensors.includeAllNetExternalWrenchesAsSensors          = false;
berdyOptionsRedundantSensors.includeAllNetExternalWrenchesAsDynamicVariables = false;
berdyOptionsRedundantSensors.includeAllJointAccelerationsAsSensors           = true;
berdyOptionsRedundantSensors.includeAllJointTorquesAsSensors                 = true;
% berdyOptionsRedundantSensors.jointOnWhichTheInternalWrenchIsMeasured.push_back('joint_1_2');
berdy_red = iDynTree.BerdyHelper();
berdy_red.init(mdlLoader.model(),mdlLoader.sensors(),berdyOptionsRedundantSensors);


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
buffers.did = iDynTree.VectorDynSize();
berdy.resizeAndZeroBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
berdy.getBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
D  = buffers.Did.toMatlab();
bD = buffers.bDid.toMatlab();
Y  = buffers.Yid.toMatlab();
bY = buffers.bYid.toMatlab();

nrOfBerdyDynEquations = size(bD,1);
[nrOfBerdySensors,nrOfBerdyDynVariables] = size(Y);
buffers.did.resize(nrOfBerdyDynVariables);

% Generate training trajectory 
% parameters 
dofs = berdy.model().getNrOfDOFs();
params.freq = 2.0; 
params.A    = 1.0;
params.w    = repmat(2*pi*params.freq,dofs,1);
genRelNoise = 0.2;
params.ddqRelNoise = genRelNoise;
params.trqsRelNoise = genRelNoise;
params.accRelNoise = genRelNoise;
params.ftRelNoise = genRelNoise;
params.stdDevDynEq  = 1e-3;
params.stdDevDynVariablesPrior = 1e3;

dofs = berdy.model().getNrOfDOFs();
% dataset.t = 0:0.02:2;
% nrOfSamples = length(dataset.t);
nrOfSamples = 15;
dataset.t   = 1:nrOfSamples;

% dataset.q = sin(params.w*dataset.t);
% dataset.dq = diag(params.w)*cos(params.w*dataset.t);
% dataset.ddq = -diag(params.w.*params.w)*sin(params.w*dataset.t);

dataset.q = rand(dofs,nrOfSamples);
dataset.dq = rand(dofs,nrOfSamples);
dataset.ddq = rand(dofs,nrOfSamples);

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

% Compute datset.acc and dataset.ft using BerdyHelpers 
dataset.acc  = zeros(3,nrOfSamples);
dataset.ft   = zeros(6,nrOfSamples);


for i = 1:nrOfSamples
    buffers.jntPos.fromMatlab(dataset.q(:,i));
    buffers.jntVel.fromMatlab(dataset.dq(:,i));
    buffers.jntAcc.fromMatlab(dataset.ddq(:,i));
    extWrenches = iDynTree.LinkWrenches();
    extWrenches.resize(berdy_red.model());
    
    berdy_red.updateKinematicsFromTraversalFixedBase(buffers.jntPos,buffers.jntVel,buffers.grav);
    berdy_red.getBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
    berdy_red.serializeDynamicVariablesComputedFromFixedBaseRNEA(buffers.jntAcc,extWrenches,buffers.did);

    Y = buffers.Yid.toMatlab();
    bY = buffers.bYid.toMatlab();
    d  = buffers.did.toMatlab();

    expectedMeas = Y*d+bY;
    dataset.acc(:,i) = expectedMeas(1:3);
    dataset.ft(:,i)  = expectedMeas(end-5,end);
    assert(length(expectedMeas) == 3+2*dofs);
end

% Generate simulated measurements, for now just using jnt torques and joint
% accelerations
measurements.ddq = zeros(size(dataset.ddq));
measurements.trqs = zeros(size(dataset.trqs));

params.stdDev_ddq = params.ddqRelNoise*std(dataset.ddq')';
params.stdDev_trqs = params.trqsRelNoise*std(dataset.trqs')';
params.stdDev_acc = params.accRelNoise*std(dataset.acc')';
params.stdDev_ft = params.ftRelNoise*std(dataset.ft')';
for i = 1:nrOfSamples
    measurements.ddq(:,i) = dataset.ddq(:,i) + normrnd(0,params.stdDev_ddq,dofs,1);
    measurements.trqs(:,i) = dataset.trqs(:,i) + normrnd(0,params.stdDev_trqs,dofs,1);
    measurements.acc(:,i) = dataset.acc(:,i) + normrnd(0,params.stdDev_acc,length(params.stdDev_acc),1);
    measurements.ft(:,i) = dataset.ft(:,i) + normrnd(0,params.stdDev_ft,length(params.stdDev_ft),1);
end

paramsStr = sprintf('ddq %.2e trqs %.2e dynEq %.2e dynVar %.2e nrOfSamples %d',params.stdDev_ddq,params.stdDev_trqs,params.stdDevDynEq,params.stdDevDynVariablesPrior,nrOfSamples);

%% Compute identifiable parameters from the specified dataset (currently hardcode to always give just the more relevant parameter)
identifibleParamsMatrix = computeJointTorquesExperimentalIdentifiableParams(dynComp,buffers,params,dataset,measurements);
[nrOfIdentifiableParams,nrOfTotalParams] = size(identifibleParamsMatrix);
dynComp.getModelDynamicsParameters(buffers.fullCadParams);
groundTruthValues = identifibleParamsMatrix*buffers.fullCadParams.toMatlab();

% For simplifyng the log data exploration, we make sure that the selected
% parameters is positive 
if( groundTruthValues(1) < 0.0 )
    identifibleParamsMatrix = -identifibleParamsMatrix;
end
groundTruthValues = identifibleParamsMatrix*buffers.fullCadParams.toMatlab();

%% Load covariances matrices 
covs.cov_e_given_d = params.stdDevDynEq.*params.stdDevDynEq*eye(nrOfBerdyDynEquations,nrOfBerdyDynEquations);
covs.cov_d         = params.stdDevDynVariablesPrior*params.stdDevDynVariablesPrior*eye(nrOfBerdyDynVariables,nrOfBerdyDynVariables);
% Sensors should be 2*dofs : first the accelerations then the torques 
% This should be modified to have a nice introspection mechanism, but for
% now this is ok 
assert(nrOfBerdySensors == 2*berdy.model().getNrOfDOFs());
nrOfDOFs = berdy.model().getNrOfDOFs();
cov_ddq = diag(params.stdDev_ddq.*params.stdDev_ddq);
cov_trqs = diag(params.stdDev_ddq.*params.stdDev_trqs);
covs.cov_y_given_d   = blkdiag(cov_ddq,cov_trqs);
measurements.y = [measurements.ddq;measurements.trqs];

%% Load covariances matrices for redundant sensors
% Sensors serialization is accelerometers, joint_accs, joint_trqs, base f/t  (again, this should 
% be exposed with a nice introspection mechanism)
covs_red = covs;
measurements_red = measurements;
cov_acc = diag(params.stdDev_acc.*params.stdDev_acc);
cov_ft = diag(params.stdDev_ft.*params.stdDev_ft);
covs_red.cov_y_given_d   = blkdiag(cov_acc,cov_ddq,cov_trqs);
measurements_red.y = [measurements.acc;measurements.ddq;measurements.trqs];

%% Let's define the loglikeliwood
llOnSomeSamples = @(inertialParams,nrOfSamples) getLogLikeLihoodOfInertialParametersUsingALimitedSetOfSamples(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,inertialParams,nrOfSamples);
llOnAllSamples = @(inertialParams) getLogLikeLihoodOfInertialParameters(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,inertialParams);
minusllOnAllSamples = @(inertialParams) -llOnAllSamples(inertialParams);

llOnSomeSamplesRed = @(inertialParams,nrOfSamples) getLogLikeLihoodOfInertialParametersUsingALimitedSetOfSamples(berdy_red,buffers,covs_red,dataset,measurements_red,identifibleParamsMatrix,inertialParams,nrOfSamples);
llOnAllSamplesRed = @(inertialParams) getLogLikeLihoodOfInertialParameters(berdy_red,buffers,covs_red,dataset,measurements_red,identifibleParamsMatrix,inertialParams);
minusllOnAllSamplesRed = @(inertialParams) -llOnAllSamplesRed(inertialParams);

% Get ll of actual  data 
llGroundTruth = getLogLikeLihoodOfInertialParameters(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,groundTruthValues);
% 
% % Let's explore the face of the ll function (not in the case of multiple links)
% inParBoundsMin = abs(groundTruthValues(2))*-5.0;
% inParBoundsMax = abs(groundTruthValues(2))*5.0;
% inPar = inParBoundsMin:0.1:inParBoundsMax;
% llExploration = zeros(size(inPar));
% for sampleIdx = 1:size(inPar,2);
%     fprintf('Computing ll for inertial parameters %f (%d out of %d)\n',inPar(sampleIdx),sampleIdx,size(inPar,2))
%     testedInParams = groundTruthValues;
%     testedInParams(2) = inPar(sampleIdx);
%     llExploration(sampleIdx) = llOnAllSamples(testedInParams);
% end
% 
% figure;
% plot(inPar,(llExploration),'b.');
% hold on;
% dim = [.2 .5 .3 .3];
% annotation('textbox',dim,'String',paramsStr,'FitBoxToText','on');
% title('LogLikelihood function');

%% After we saw the shape of the LL, lets 
identification_enabled = true;

if( identification_enabled )
    fprintf('Identification enabled, tryng to identify the inertial parameters');
    %% Estimate inertial parameters, return cell array of estimated parameters for number of used samples 
    classicalEstimationResults = getClassicalEstimationResults(dynComp,buffers,dataset,measurements,identifibleParamsMatrix);

    %% Estimate inertial parameters using Berdy maximization of LL 
    fprintf('Maximizing ll on the complete set of partial datasets\n');
    tic;
    berdyEstimationResults   = getBerdyEstimationResults(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,groundTruthValues);
    timeSpent = toc;
    fprintf('Seconds spent in optimizing complete set of datasets LL : %f\n',timeSpent);
    
    %% Estimate inertial parameters using Berdy maximization of LL 
    % fprintf('Maximizing ll on the complete set of partial datasets\n');
    % tic;
    % berdyRedEstimationResults   = getBerdyEstimationResults(berdy_red,buffers,covs_red,dataset,measurements_red,identifibleParamsMatrix,groundTruthValues);
    % timeSpent = toc;
    % fprintf('Seconds spent in optimizing complete set of datasets LL : %f\n',timeSpent);

    %% Plot relative error in estimation results 
    classicalEstimationErrors = (classicalEstimationResults-repmat(groundTruthValues,1,nrOfSamples));
    classicalEstimationRelativeErrorsNorms = sqrt(sum(abs(classicalEstimationErrors).^2,1))./norm(groundTruthValues);

    berdyEstimationErrors = (berdyEstimationResults-repmat(groundTruthValues,1,nrOfSamples));
    berdyEstimationRelativeErrorsNorms = sqrt(sum(abs(berdyEstimationErrors).^2,1))./norm(groundTruthValues);
    
    % berdyRedEstimationErrors = (berdyRedEstimationResults-repmat(groundTruthValues,1,nrOfSamples));
    % berdyRedEstimationRelativeErrorsNorms = sqrt(sum(abs(berdyRedEstimationErrors).^2,1))./norm(groundTruthValues);

    figure;
    lineWidth = 4;
    plot(dataset.t(1:end),classicalEstimationRelativeErrorsNorms,'r--','LineWidth',lineWidth);
    hold on;
    plot(dataset.t(1:end),berdyEstimationRelativeErrorsNorms,'b-.','LineWidth',lineWidth);
    % plot(dataset.t(1:end),berdyRedEstimationRelativeErrorsNorms,'b=.-');
    plot(dataset.t(1:end),zeros(size(dataset.t(1:end))),'g','LineWidth',lineWidth);
    %dim = [.2 .5 .3 .3];
    %annotation('textbox',dim,'String',paramsStr,'FitBoxToText','on');
    title('4 links inertial parameters identification')
    xlabel('# of considered samples');
    ylabel('Relative norm error');
    legend('Least Squares','Likelihood Maximization','Ground Truth');
    set(gca,'FontSize',25);
    
end

