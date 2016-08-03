function classicalEstimationResults = getLogLikeLihoodOfInertialParameters(dynComp, buffers, dataset, measurements, identifibleParamsMatrix, baseParams, berdyOptions)
% Given a set of inertial parameters and a set of sensors return its log likeliwood over all the dataset  
% check https://www.sharelatex.com/project/57a090a9751e872d5f995acb for the
% formulas
[dofs,nrOfSamples] = size(dataset.dq);
[nrOfBaseParameters,nrOfTotalParameters] = size(identifibleParamsMatrix);

% Create a new berdy with the specified options and the passed
% inertialParams 
newInertialParams = integrateIdentifiedParams(baseParams,cadParams,identifibleParamsMatrix);



classicalEstimationResults = zeros(nrOfBaseParameters,nrOfSamples);

for i = 1:nrOfSamples
    trqRegr = getTrqsBaseRegressor(dynComp, buffers, dataset.q(:,i), dataset.dq(:,i), measurements.ddq(:,i), identifibleParamsMatrix);
    accRegrs = accRegrs + trqRegr'*trqRegr;
    accKnownTerms = accKnownTerms + trqRegr'*measurements.trqs(:,i);
    
    classicalEstimationResults(:,i) = pinv(accRegrs)*accKnownTerms;
end

end

