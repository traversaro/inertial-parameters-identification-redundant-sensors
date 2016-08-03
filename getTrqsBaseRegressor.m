function [ trqDynBaseRegr ] = getTrqsBaseRegressor(dynComp, buffers, q, dq, ddq , identifibleParamsMatrix )
    trqDynBaseRegr = getTrqsRegressor(dynComp,buffers,q,dq,ddq)*identifibleParamsMatrix';
end

