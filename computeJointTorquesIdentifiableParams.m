function [ idParamsMatrix ] = computeJointTorquesIdentifiableParams( dynComp, state, params )
%computeIdentifiableParams Compute inertial parameters identifiable from
%fixed base torque measurements

dofs = dynComp.getNrOfDegreesOfFreedom();
grav = [0.0,0.0,-10.0,0.0,0.0,0.0];

qId = iDynTree.VectorDynSize();
dqId = iDynTree.VectorDynSize();
ddqId = iDynTree.VectorDynSize();
gravId = iDynTree.SpatialAcc();
gravId.fromMatlab(grav);
fullDynRegrId = iDynTree.MatrixDynSize();
fullDynRegrId.zero();
nrOfInertialParams = 10*dynComp.getNrOfLinks();
acc = zeros(nrOfInertialParams,nrOfInertialParams);

for i = 1:100
    q = rand(dofs,1);
    dq = rand(dofs,1);
    ddq = rand(dofs,1);
    
    trqDynRegr = getTrqsRegressor( dynComp, state, q, dq, ddq);
        
    % Accumulate
    acc = acc + trqDynRegr'*trqDynRegr;
end

[U,S,V] = svd(acc);
idParamsMatrix = (V(:,1:rank(acc)))';

end

