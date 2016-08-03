function [ trqDynRegr ] = getTrqsRegressor(dynComp, buffers, q, dq, ddq )
    buffers.jntPos.fromMatlab(q);
    buffers.jntVel.fromMatlab(dq);
    buffers.jntAcc.fromMatlab(ddq);
    
    dynComp.setRobotState(buffers.jntPos,buffers.jntVel,buffers.jntAcc,buffers.gravSpatialAcc);
    dynComp.getDynamicsRegressor(buffers.fullDynRegrId);
    
    % Extract the part relative to torques 
    fullDynRegr = buffers.fullDynRegrId.toMatlab();
    trqDynRegr = fullDynRegr(7:end,:);
end

