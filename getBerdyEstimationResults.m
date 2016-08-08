function [berdyEstimationResults] = getBerdyEstimationResults(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,groundTruthValues)
[dofs,nrOfSamples] = size(dataset.dq);
[nrOfBaseParameters,nrOfTotalParameters] = size(identifibleParamsMatrix);

berdyEstimationResults = zeros(nrOfBaseParameters,nrOfSamples);

for i = 1:nrOfSamples
    fprintf('Running Berdy ll maximization for samples from 1 to %d (out of %d)\n',i,nrOfSamples);
    % Create ll including the samples from 1 to 1 
    minusll = @(pi) -getLogLikeLihoodOfInertialParametersUsingALimitedSetOfSamples(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,pi,i);
    
    % We assume that we have meaningful a-priori CAD parameters to start the optimization process 
    startingPoint = 0.95.*groundTruthValues;
    
    % Find the minimux of minusll (i.e. the max of ll)
    options = optimset('MaxFunEvals', 1000000);
    options = optimset('MaxIter', 200000);
    berdyEstimationResults(:,i) = fminsearch(minusll,startingPoint,options);
end

end

