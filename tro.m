% Load model
mdlLoader = iDynTree.ModelLoader();
mdlLoader.loadModelFromFile('threeLinks.urdf');

% Create berdyHelper 
berdyOptions = iDynTree.BerdyOptions();
berdy = iDynTree.BerdyHelper();
berdy.init(mdlLoader.model(),mdlLoader.sensors(),berdyOptions);

% Create state 
state.jntPos = iDynTree.JointPosDoubleArray(berdy.model());
state.jntVel = iDynTree.JointDOFsDoubleArray(berdy.model());
state.jntPos.zero();
state.jntVel.zero();

% Gravity 
state.grav = iDynTree.Vector3();
state.grav.fromMatlab([0.0 0.0 -10.0]);

% Update kinematics 
berdy.updateKinematicsFromTraversalFixedBase(state.jntPos,state.jntVel,state.grav);

% Get matrices 
Yid = iDynTree.MatrixDynSize();
Did = iDynTree.MatrixDynSize();
bYid = iDynTree.VectorDynSize();
bDid = iDynTree.VectorDynSize();
berdy.resizeAndZeroBerdyMatrices(Did,bDid,Yid,bYid);
berdy.getBerdyMatrices(Did,bDid,Yid,bYid);
D  = Did.toMatlab();
bD = bDid.toMatlab();
Y  = Yid.toMatlab();
bY = bYid.toMatlab();
