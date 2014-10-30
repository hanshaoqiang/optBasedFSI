%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
%
% Task : Main script that executes the FSI simulation.
%

%% Preamble
clc;
clear;
close all;

%% Includes
addpath('../mesh');
addpath('../mesh/Mesh_Turek_RigidCylinder');
casename = 'Turek_FSI_BM_RigidCyl';
% addpath('../mesh/Mesh_Turek_RigidCylinder');
addpath('../solver');
addpath('../interface');
addpath('../graphics');
addpath('../graphics/Output_vtk');
addpath('../fsi_coupling/');
addpath('../empire/');


%% For FLUID Solver
% Add all functions related to general graphics of post processing
addpath('../graphics/fluid_Graphics');

% Path for Mesh movement solver
addpath('../fluidSolver/possionEquationSolver');
% Adding the steady state solvers
addpath('../steadyStateSolvers/fluid');


% Add all functions related to the Computational Fluid
% Dynamics problems
addpath('../fluidSolver/femForCFD/stiffnessMatrices',...
    '../fluidSolver/femForCFD/solvers',...
    '../fluidSolver/femForCFD/neumannBoundaryConditions',...
    '../fluidSolver/femForCFD/graphics',...
    '../fluidSolver/femForCFD/postprocessing',...
    '../fluidSolver/femForCFD/massMatrices',...
    '../fluidSolver/femForCFD/inhomogeneousDirichletBoundaryConditions',...
    '../fluidSolver/femForCFD/initialConditions',...
    '../fluidSolver/mathematics');



%% For STRUCTURAL Solver
% Add all functions related to general graphics of post processing
addpath('../graphics/structure_Graphics');


% Add all functions related to the Computational Fluid
% Dynamics problems
addpath('../structuralSolver/basisFunctions',...
    '../structuralSolver/loads',...
    '../structuralSolver/postprocessing',...
    '../structuralSolver/solvers',...
    '../structuralSolver/stiffnessMatrices',...
    '../structuralSolver/supports');


% Adding the steady state solvers
addpath('../steadyStateSolvers/structure');


%% For MESH MOTION Solver
% Add all functions related to general graphics of post processing
addpath('../graphics/structure_Graphics');


% Add all functions related to the Computational Fluid
% Dynamics problems
addpath('../meshSolver/basisFunctions',...
    '../meshSolver/loads',...
    '../meshSolver/postprocessing',...
    '../meshSolver/solvers',...
    '../meshSolver/stiffnessMatrices',...
    '../meshSolver/supports');



%% Mean velocity U bar of the inlet parabolic flow
Ubar = 2.0;


%% Physical properties of Fluids

% Kinematic viscosity
fphysics.nue = 1e-3;

% Density of Fluid
fphysics.rho = 1e3;

% Source vector
amplification = 0;
b = amplification*[1 0]';

% Definition of body Force on the fluid
fphysics.bodyForce.magnitude = 0;
fphysics.bodyForce.direction = [1;0;0];


% Integration scheme
% type = 0 : default FGI integration element-wise
% type = 1 : manual choice of the number of Gauss points
fphysics.int.type = 0;
% Integrating the stiffness
fphysics.int.nGauss = 1;


% On the stabilization
% stabilization.type = 'noStabilization';
fphysics.stabilization.type = 'variationalMultiscale';
fphysics.stabilization.parametersType = 'automatic';
fphysics.stabilization.CI = 4;
fphysics.stabilization.Ct = 1;
fphysics.stabilization.scalingTauM = 1;
fphysics.stabilization.scalingTauC = 1;


% On the iterations for the complete nonlinear system
fphysics.nonLinearScheme.method = 'Newton';
fphysics.nonLinearScheme.eps = 1e-6;
fphysics.nonLinearScheme.maxIter = 15;


%% Physical properties of Structure

% Young's modulus
sphysics.E = 0.50*1e6;

% Poisson ratio
sphysics.nu = 0.4;

% Density of the material
sphysics.rho = 1e3;

% Thickness of the plate
sphysics.t = 1;

% Number of Integration points
sphysics.int.nGauss = 3;

% Body force direction
sphysics.bodyForce = [0,-0];

%% GUI

% Analysis type
sphysics.dimension = '2d';
sphysics.dofs = 'displacements';
sphysics.physics = 'plainStrain';
sphysics.type = 'linear';


%% The start and the end time of the simulation
transient.Tstart = 0;
transient.Tend = 10;

% The number of time steps
transient.nT = 2000;

% The time step
transient.dt = (transient.Tend - transient.Tstart)/(transient.nT);

% Parameters to ensure stability
% beta > 0.25
transient.beta = 0.5;

% alphaBeta <= .5
transient.alphaBeta = -0.1;

% alphaBeta + gamma >= .25
transient.gamma = .5 - transient.alphaBeta;

% The frequency at which the update step between structure and fluid takes
% places
transient.updateFreq = 1;
transient.noFluidItrFSI = 1;
transient.noFluidItrFSI_inFPI = 1;

% Residual for stopping the fixed point iterations of FSI solver
transient.fsiProp.fsiResidual = 1e-9;
% Maximum number of fixed point iterations
transient.fsiProp.maxFpIterations = 35;
% Type of residual for constant relaxation.
%  1 --- based on force on fsi interface
%  2 --- based on displacement of fsi interface
transient.fsiProp.residualType = 2;

% The entitiy on which relaxation is done during the Aitken relaxation
% 1 -- on Forces
% 2 -- on Structure displacements
transient.fsiProp.aitkenResidualType = 2;

% Relaxation factor for the Fixed point iterations
transient.fsiProp.omega = 0.4;

%% Parameters for ale mesh movement
physics_ale = sphysics;
physics_ale.rho = 1;

%% Obtaining the mesh and boundary conditions for both FLUID and STRUCTURE from GID 
[mesh_s, bc_s, mesh_f, bc_f, mesh_ale, bc_ale] = readMeshFromGID(casename,Ubar);

%% Initializing the EMPIRE environment
EMPIRE_CORE_BASE_DIR=getenv('EMPIRE_CORE_BASE_DIR');
mexFuncPath=strcat(EMPIRE_CORE_BASE_DIR, '/interface/matlab/binICC');
path(path, mexFuncPath);
EMPIRE_API_Connect('matlabClient.xml');
display(mexFuncPath);
% format specification in EMPIRE/EMPIRE-Testing/benchmark/FieldExamples/matlabExample/matlabClient.m
empireMesh = obtainMeshToSendToEMPIRE(mesh_s); % Obtain the interface mesh to send to EMPIRE Format should be followed.
EMPIRE_API_sendMesh('',empireMesh.numNodes, empireMesh.numElems, empireMesh.nodes,...
                     empireMesh.nodeIDs, empireMesh.numNodesPerElem, empireMesh.elems); % Send mesh to EMPIRE

%% Solving the FSI problem by a call to the solver
solve_FSIProblemBossakScheme2D(mesh_f, bc_f, fphysics, mesh_s, bc_s, sphysics, mesh_ale, bc_ale, physics_ale, transient);


% End of the script
