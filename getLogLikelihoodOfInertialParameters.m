function [logLikeLihood, likelihood]  = getLogLikeLihoodOfInertialParameters(berdy, buffers, dataset, measurements, identifibleParamsMatrix, baseParams)
% Given a set of inertial parameters and a set of sensors return its log likeliwood over all the dataset  
% check https://www.sharelatex.com/project/57a090a9751e872d5f995acb for the
% formulas
[dofs,nrOfSamples] = size(dataset.dq);
[nrOfBaseParameters,nrOfTotalParameters] = size(identifibleParamsMatrix);

% Save the inertial parameters in input 
berdy.model().getInertialParameters(buffers.fullCadParams);

% Compute a suitable set of inertial parameters 
newInertialParams = integrateIdentifiedParams(baseParams,buffers.fullCadParams.toMatlab(),identifibleParamsMatrix);
buffers.fullEstimatedParams.fromMatlab(newInertialParams);

% Load this new inertial parameters in Berdy's model
berdy.model().updateInertialParameters(buffers.fullEstimatedParams);

accLogLikeLihood = 0;

for i = 1:nrOfSamples
    accLogLikeLihood = accLogLikeLihood + getLLOfSingleSample(berdy, buffers, dataset, measurements);
end

% Restore the original inertial parameters in berdy 
berdy.model().updateInertialParameters(buffers.fullCadParams);

logLikeLihood = accLogLikeLihood;
likelihood = exp(logLikeLihood);

end

