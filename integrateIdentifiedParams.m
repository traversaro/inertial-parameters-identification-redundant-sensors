function newInertialParams = integrateIdentifiedParams(baseParams,cadParams,identifibleParamsMatrix)
% Integrate the identified parameters in available cadParams
% I.e. use the identified baseParams on the identifiable subspace, while 
% on the non-identifiable subspace use the projection of the CAD params 
    fullParametersSize = size(cadParams,1);
    newInertialParams = identifibleParamsMatrix'*baseParams + (eye(fullParametersSize,fullParametersSize)-identifibleParamsMatrix'*identifibleParamsMatrix)*cadParams;