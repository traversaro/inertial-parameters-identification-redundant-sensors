function [ idParamsMatrix ] = computeJointTorquesExperimentalIdentifiableParams( dynComp, buffers, params , dataset, measurements)
%computeIdentifiableParams Compute inertial parameters identifiable from
%fixed base torque measurements

nrOfInertialParams = 10*dynComp.getNrOfLinks();
[dofs,nrOfSamples] = size(dataset.dq);

acc = zeros(nrOfInertialParams,nrOfInertialParams);

for i = 1:nrOfSamples
    q = dataset.q(:,i);
    dq = dataset.dq(:,i);
    ddq = measurements.ddq(:,i);
    
    trqDynRegr = getTrqsRegressor( dynComp, buffers, q, dq, ddq);
        
    % Accumulate
    acc = acc + trqDynRegr'*trqDynRegr;
end

[U,S,V] = svd(acc);
idParamsMatrix = (V(:,1))';

end

