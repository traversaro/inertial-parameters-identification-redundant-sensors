function classicalEstimationResults = getClassicalEstimationResults(dynComp,buffers,dataset,measurements,identifibleParamsMatrix)
[dofs,nrOfSamples] = size(dataset.dq);
[nrOfBaseParameters,nrOfTotalParameters] = size(identifibleParamsMatrix);

accRegrs = zeros(nrOfBaseParameters,nrOfBaseParameters);
accKnownTerms = zeros(nrOfBaseParameters,1);

classicalEstimationResults = zeros(nrOfBaseParameters,nrOfSamples);

for i = 1:nrOfSamples
    trqRegr = getTrqsBaseRegressor(dynComp, buffers, dataset.q(:,i), dataset.dq(:,i), measurements.ddq(:,i), identifibleParamsMatrix);
    accRegrs = accRegrs + trqRegr'*trqRegr;
    accKnownTerms = accKnownTerms + trqRegr'*measurements.trqs(:,i);
    
    classicalEstimationResults(:,i) = pinv(accRegrs)*accKnownTerms;
end

end

