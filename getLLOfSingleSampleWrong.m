function [ llSingleSample ] = getLLOfSingleSample(berdy, buffers, covs, dataset, measurements , sampleIdx  )
    buffers.jntPos.fromMatlab(dataset.q(:,sampleIdx));
    buffers.jntVel.fromMatlab(dataset.dq(:,sampleIdx));
    berdy.updateKinematicsFromTraversalFixedBase(buffers.jntPos,buffers.jntVel,buffers.grav);
    
    berdy.getBerdyMatrices(buffers.Did,buffers.bDid,buffers.Yid,buffers.bYid);
    D  = buffers.Did.toMatlab();
    bD = buffers.bDid.toMatlab();
    Y  = buffers.Yid.toMatlab();
    bY = buffers.bYid.toMatlab();
    
    % Let's roll with some computations, check https://www.sharelatex.com/project/57a090a9751e872d5f995acb
    Sigma_D = covs.cov_e_given_d + D*covs.cov_d*D';
    invSigma_D = inv(Sigma_D);
    % Sigma_overline_d = covs.cov_d+covs.cov_d*D'*invSigma_D*D*covs.cov_d;
    Sigma_overline_d   = inv(D'*inv(covs.cov_e_given_d)*D+inv(covs.cov_d));
    Sigma_y_ll = covs.cov_y_given_d + Y*Sigma_overline_d*Y';
    mu_y_ll    = bY + Y*(covs.cov_d*D'*invSigma_D*bD);
    detSigma_y_ll = det(Sigma_y_ll);
    nrOfMeasurements = size(Y,1);
    llSingleSample = -(measurements.y(:,sampleIdx)-mu_y_ll)'*(Sigma_y_ll\(measurements.y(:,sampleIdx)-mu_y_ll))-0.5*log(detSigma_y_ll);%-0.5*nrOfMeasurements*log(2*pi);
end

