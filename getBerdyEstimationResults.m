function [berdyEstimationResults] = getBerdyEstimationResults(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix)
[dofs,nrOfSamples] = size(dataset.dq);
[nrOfBaseParameters,nrOfTotalParameters] = size(identifibleParamsMatrix);

berdyEstimationResults = zeros(nrOfBaseParameters,nrOfSamples);

for i = 1:nrOfSamples
    fprintf('Running berding ll maximization for samples from 1 to %d\n',i);
    % Create ll including the samples from 1 to 1 
    minusll = @(pi) -getLogLikeLihoodOfInertialParametersUsingALimitedSetOfSamples(berdy,buffers,covs,dataset,measurements,identifibleParamsMatrix,pi,i);
    
    if( i == 1 )
        startingPoint = 0.2*ones(nrOfBaseParameters,1);
    end
    
    % Find the minimux of minusll (i.e. the max of ll)
    berdyEstimationResults(:,i) = fminsearch(minusll,startingPoint);
end

end

