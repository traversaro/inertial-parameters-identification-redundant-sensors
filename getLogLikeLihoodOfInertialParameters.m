function [logLikeLihood, likelihood]  = getLogLikeLihoodOfInertialParameters(berdy, buffers, covs, dataset, measurements, identifibleParamsMatrix, baseParams)
    [~,nrOfSamples] = size(dataset.dq);
    [logLikeLihood, likelihood]  = getLogLikeLihoodOfInertialParametersUsingALimitedSetOfSamples(berdy, buffers, covs, dataset, measurements, identifibleParamsMatrix, baseParams, nrOfSamples);
end

