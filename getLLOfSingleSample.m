function [ llSingleSample ] = getLLOfSingleSample(berdy, buffers, covs, dataset, measurements , sampleIdx  )
    buffers.jntPos.fromMatlab(dataset.q(:,sampleIdx));
    buffers.jntVel.fromMatlab(dataset.dq(:,sampleIdx));
    berdy.updateKinematicsFromTraversalFixedBase(buffers.jntPos,buffers.jntVel,buffers.grav);
    
    berdy.getBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
    D  = buffers.Did.toMatlab();
    bD = buffers.bDid.toMatlab();
    Y  = buffers.Yid.toMatlab();
    bY = buffers.bYid.toMatlab();
    
    % Let's roll with some computations, check http://wiki.codyco.eu/doku.php?id=analytic_solution_of_the_maximum-a-posteriori_dynamics
    % overlineSigma_D = inv(inv(covs.cov_d)+D'*inv(covs.cov_e_given_d)*D);
    overlineSigma_D = covs.cov_d - covs.cov_d*D'*inv(D*covs.cov_d*D'+covs.cov_e_given_d)*D*covs.cov_d;
    overlineMu_D = -overlineSigma_D*(D'*inv(covs.cov_e_given_d)*bD);
    
    Sigma_y_ll = covs.cov_y_given_d + Y*overlineSigma_D*Y';
    mu_y_ll    = bY + Y*overlineMu_D;
    detSigma_y_ll = det(Sigma_y_ll);
    nrOfMeasurements = size(Y,1);
    llSingleSample = -(measurements.y(:,sampleIdx)-mu_y_ll)'*(Sigma_y_ll\(measurements.y(:,sampleIdx)-mu_y_ll))-0.5*log(detSigma_y_ll);%-0.5*nrOfMeasurements*log(2*pi);
end

