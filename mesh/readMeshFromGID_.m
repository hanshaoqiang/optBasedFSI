%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [mesh_s, bc_s, mesh_f, bc_f, mesh_ale, bc_ale] = readMeshFromGID
%% This function reads the output files form the GID output and formulates
% the necessary arrays and element tables for both fluid and structural
% solvers.
%
% The mesh files are by default assumed to be in the path
% /mesh/GID_output
%
%
%   OUTPUT
%
%       mesh_s      : mesh with entities for structrual solver
%       mesh_f      : mesh for entities for fluid solver
%       bc_s        : boundary conditions for structural problem
%       bc_f        : boundary conditions for fluid solver
%       mesh_ale    : mesh required for ALE formulation essentially the
%                     same as fluid mesh but with different boundary
%                     conditions
%       bc_ale      : boundary conditions for ALE Mesh Solver
%


%% Loading all the output files of the GID mesh
load('fluidNodes.txt');
load('structureNodes.txt');
load('fluidElements.txt');
load('structureElements.txt');

load('fluidWallNodes.txt');
load('fluidInletNodes.txt');
load('fluidOutletNodes.txt');

load('fsiInterfaceNodes.txt');

load('structureFixedNodes.txt');
load('structureForceNodes.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% FLUID %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Formulating entities for fluid mesh
mesh_f.nodes = fluidNodes(:,2:4);
mesh_f.nodesGlobalNumbers  = fluidNodes(:,1);
mesh_f.elements = fluidElements(:,2:4);
mesh_f.inletNodes = fluidInletNodes(:,1);
mesh_f.wallNodes = fluidWallNodes(:,1);
mesh_f.outletNodes = fluidOutletNodes(:,1);
mesh_f.fsiNodes = fsiInterfaceNodes(:,1);
mesh_f.fsiNodesGlobal = mesh_f.fsiNodes;
% Initial mesh velocity for ALE purposes
mesh_f.velocity = zeros(3*length(mesh_f.nodes),1);

% Making the above read global (including fluid and structure) numbers in
% to local means nodes start at 1
% Loop over all the elements to find what are the local numbers of its
% nodes
for i = 1:length(mesh_f.elements)
    node1_global = mesh_f.elements(i,1);
    node2_global = mesh_f.elements(i,2);
    node3_global = mesh_f.elements(i,3);
    
    node1_local = find(node1_global == mesh_f.nodesGlobalNumbers);
    node2_local = find(node2_global == mesh_f.nodesGlobalNumbers);
    node3_local = find(node3_global == mesh_f.nodesGlobalNumbers);
    
    mesh_f.elements(i,1) = node1_local;
    mesh_f.elements(i,2) = node2_local;
    mesh_f.elements(i,3) = node3_local;
end

% Doing the same thing for inletNodes
for i = 1:length(mesh_f.inletNodes)
    node_global = mesh_f.inletNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.inletNodes(i) = node_local;
end

% For wallNodes
for i = 1:length(mesh_f.wallNodes)
    node_global = mesh_f.wallNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.wallNodes(i) = node_local;
end

% For outletNodes
for i = 1:length(mesh_f.outletNodes)
    node_global = mesh_f.outletNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.outletNodes(i) = node_local;
end

% For fsiNodes
for i = 1:length(mesh_f.fsiNodes)
    node_global = mesh_f.fsiNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.fsiNodes(i) = node_local;
end


%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
bc_f.drichletDOF = [3*mesh_f.inletNodes-2; 3*mesh_f.inletNodes-1; 3*mesh_f.wallNodes-2; 3*mesh_f.wallNodes-1; 3*mesh_f.fsiNodes-2; 3*mesh_f.fsiNodes-1];

% Drichlet values vector
bc_f.drichletVector = zeros(3*length(mesh_f.nodes),1);
bc_f.drichletVector(3*mesh_f.inletNodes-2) = fluidInletNodes(:,2); % Inlet X velocity
bc_f.drichletVector(3*mesh_f.inletNodes-1) = fluidInletNodes(:,3); % Inlet Y velocity
% IMPORTANT :: IMPORTANT :: Keep in mind that the fsi Values should be
% updated at every time step


% Neumann DOF vector
bc_f.neumannDOF = 3*mesh_f.outletNodes;
% Neumann values vector
% For outlet condition we do have a 0 neumann condition
bc_f.neumannVector = zeros(3*length(mesh_f.nodes),1);

% Finding the DOFs which are to be computed
bc_f.freeDOF = setdiff(1:3*length(mesh_f.nodes), bc_f.drichletDOF);

% Formulating the FSI nodes on fluid side
mesh_f.fsiDOF = zeros(2*length(mesh_f.fsiNodes),1);
mesh_f.fsiDOF(1:2:end) = 3*mesh_f.fsiNodes-2;
mesh_f.fsiDOF(2:2:end) = 3*mesh_f.fsiNodes-1;

% mesh_f.fsiDOF = [3*mesh_f.fsiNodes-2; 3*mesh_f.fsiNodes-1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% ALE MESH  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basically all the nodes and elements are same as the fluid mesh.
mesh_ale    = mesh_f;
bc_ale      = bc_f;

%% Formulating the FSI nodes DOF
bc_ale.fsiDOF = [2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes];

%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
bc_ale.drichletDOF = [2*mesh_ale.inletNodes-1; 2*mesh_ale.inletNodes; 2*mesh_ale.wallNodes-1; 2*mesh_ale.wallNodes; 2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes; 2*mesh_ale.outletNodes-1; 2*mesh_ale.outletNodes];

% Drichlet values vector
bc_ale.drichletVector = zeros(2*length(mesh_ale.nodes),1);

% Finding the DOFs which are to be computed
bc_ale.freeDOF = setdiff(1:2*length(mesh_ale.nodes), bc_ale.drichletDOF);

% Formulating the FSI nodes on fluid side
mesh_ale.fsiDOF = zeros(2*length(mesh_ale.fsiNodes),1);
mesh_ale.fsiDOF(1:2:end) = 2*mesh_ale.fsiNodes-1;
mesh_ale.fsiDOF(2:2:end) = 2*mesh_ale.fsiNodes;
% mesh_ale.fsiDOF = [2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes];

%% ::: IMPORTANT :: IMPORTANT ::
% The FSI interface nodes Drichlet vector should always be updated in every
% time step



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %% Formulating entities for Structure mesh
mesh_s.nodes = structureNodes(:,2:4);
mesh_s.nodesGlobalNumbers  = structureNodes(:,1);
mesh_s.elements = structureElements(:,2:4);
mesh_s.fixedNodes = structureFixedNodes(:,1);
mesh_s.forceNodes = structureForceNodes(:,1);
mesh_s.fsiNodes = fsiInterfaceNodes(:,1);
mesh_s.fsiNodesGlobal = mesh_s.fsiNodes;


% Making the above read global (including fluid and structure) numbers in
% to local means nodes start at 1
% Loop over all the elements to find what are the local numbers of its
% nodes
for i = 1:length(mesh_s.elements)
    node1_global = mesh_s.elements(i,1);
    node2_global = mesh_s.elements(i,2);
    node3_global = mesh_s.elements(i,3);
    
    node1_local = find(node1_global == mesh_s.nodesGlobalNumbers);
    node2_local = find(node2_global == mesh_s.nodesGlobalNumbers);
    node3_local = find(node3_global == mesh_s.nodesGlobalNumbers);
    
    mesh_s.elements(i,1) = node1_local;
    mesh_s.elements(i,2) = node2_local;
    mesh_s.elements(i,3) = node3_local;
end

% Doing the same for fixedNodes
for i = 1:length(mesh_s.fixedNodes)
    node_global = mesh_s.fixedNodes(i);
    node_local = find(node_global == mesh_s.nodesGlobalNumbers);
    mesh_s.fixedNodes(i) = node_local;
end

% Doing the same for ForceNodes
for i = 1:length(mesh_s.forceNodes)
    node_global = mesh_s.forceNodes(i);
    node_local = find(node_global == mesh_s.nodesGlobalNumbers);
    mesh_s.forceNodes(i) = node_local;
end

% For fsiNodes
for i = 1:length(mesh_s.fsiNodes)
    node_global = mesh_s.fsiNodes(i);
    node_local = find(node_global == mesh_s.nodesGlobalNumbers);
    mesh_s.fsiNodes(i) = node_local;
end


%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
bc_s.drichletDOF = [2*mesh_s.fixedNodes-1; 2*mesh_s.fixedNodes];
% Drichlet Values vector
bc_s.drichletVector = zeros(2*length(mesh_s.nodes),1);
% % bc_s.drichletVector(2*mesh_s.fixedNodes-1) = structureFixedNodes(:,2); % Predescribed X displacement
% % bc_s.drichletVector(2*mesh_s.fixedNodes) = structureFixedNodes(:,3);   % Predescribed Y displacement
bc_s.drichletVector(2*mesh_s.fixedNodes-1) = 0.0; % Predescribed X displacement
bc_s.drichletVector(2*mesh_s.fixedNodes) = 0.0;   % Predescribed Y displacement


% Neumann DOF vector
bc_s.neumannDOF = [2*mesh_s.forceNodes-1; 2*mesh_s.forceNodes];
bc_s.neumannVector = zeros(2*length(mesh_s.nodes),1);
bc_s.neumannVector(2*mesh_s.forceNodes-1) = structureForceNodes(:,2); % X Component of force
bc_s.neumannVector(2*mesh_s.forceNodes) = structureForceNodes(:,3);   % Y Component of force

% Finding the DOFs which are to be computed
bc_s.freeDOF = setdiff(1:2*length(mesh_s.nodes), bc_s.drichletDOF);

% Formulating the FSI nodes on structure side
mesh_s.fsiDOF = zeros(2*length(mesh_s.fsiNodes),1);
mesh_s.fsiDOF(1:2:end) = 2*mesh_s.fsiNodes-1;
mesh_s.fsiDOF(2:2:end) = 2*mesh_s.fsiNodes;
% mesh_s.fsiDOF = [2*mesh_s.fsiNodes-1; 2*mesh_s.fsiNodes];

% :: IMPORTANT :: IMPORTANT :: IMPORTANT ::
% Here forces on FSI interface are not considered in the Neumann Vector.
% They should be added explicitly in the structural solver.

% Formulating the FSI interface edges on the structure side.
% mesh_s = formulateFSIedges(mesh_s);


% End of the method
end