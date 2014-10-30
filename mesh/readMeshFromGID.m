function [mesh_s, bc_s, mesh_f, bc_f, mesh_ale, bc_ale] = readMeshFromGID(casename, Ubar)
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
meshFileName = [casename '.dat'];
fstring = fileread(meshFileName); % read the mesh file as one string

%% FLUID_NODES
fblocks = regexp(fstring,'FLUID_NODES','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_f.nodes = out(:,2:4);
mesh_f.nodesGlobalNumbers = out(:,1);


%% FLUID_ELEMENTS
fblocks = regexp(fstring,'FLUID_ELEMENTS','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_f.elements = out(:,2:4);

%% WALL_NODES
fblocks = regexp(fstring,'WALL_NODES','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_f.wallNodes = out(:,1);

%% INLET_NODES
fblocks = regexp(fstring,'INLET_NODES','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_f.inletNodes = out(:,1);
% Drichlet values vector
bc_f.drichletVector = zeros(3*length(mesh_f.nodes),1);
maxVel = Ubar*norm(out(:,2),inf);


%% OUTLET_NODES
fblocks = regexp(fstring,'OUTLET_NODES','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_f.outletNodes = out(:,1);


%% STRUCTURE_NODES
fblocks = regexp(fstring,'STRUCTURE_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_s.nodes = out(:,2:4);
mesh_s.nodesGlobalNumbers = out(:,1);


%% STRUCTURE_ELEMENTS
fblocks = regexp(fstring,'STRUCTURE_ELEMENTS','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_s.elements = out(:,2:4);


%% FIXED_NODES
fblocks = regexp(fstring,'FIXED_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_s.fixedNodes = out(:,1);

%% FORCE_NODES
fblocks = regexp(fstring,'FORCE_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_s.forceNodes = out(:,1);

%% FSI_INTERFACE_NODES
fblocks = regexp(fstring,'FSI_INTERFACE_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
mesh_s.fsiNodesGlobal = out(:,1);
mesh_f.fsiNodesGlobal = out(:,1);
mesh_s.fsiNodes = out(:,1);
mesh_f.fsiNodes = out(:,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% FLUID %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Formulating entities for fluid mesh
% % % % mesh_f.nodes = fluidNodes(:,2:4);
% % % % mesh_f.nodesGlobalNumbers  = fluidNodes(:,1);
% % % % mesh_f.elements = fluidElements(:,2:4);
% % % % mesh_f.inletNodes = fluidInletNodes(:,1);
% % % % mesh_f.wallNodes = fluidWallNodes(:,1);
% % % % mesh_f.outletNodes = fluidOutletNodes(:,1);
% % % % mesh_f.fsiNodes = fsiInterfaceNodes(:,1);
% % % % mesh_f.fsiNodesGlobal = mesh_f.fsiNodes;
% Initial mesh velocity for ALE purposes
mesh_f.velocity = zeros(3*length(mesh_f.nodes),1);
mesh_f.globalWallNumbers = mesh_f.wallNodes;
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

% % % Drichlet values vector
% % bc_f.drichletVector = zeros(3*length(mesh_f.nodes),1);
% % maxVel = Ubar*norm(fluidInletNodes(:,2),inf);

y = mesh_f.nodes(mesh_f.inletNodes,2);
parabolicVelProf = 1.5*maxVel*(4.0/0.1681)*y.*(0.41-y);

bc_f.drichletVector(3*mesh_f.inletNodes-2) = parabolicVelProf; % Inlet X velocity
% bc_f.drichletVector(3*mesh_f.inletNodes-1) = fluidInletNodes(:,3); % Inlet Y velocity


% bc_f.drichletVector(3*mesh_f.inletNodes-2) = fluidInletNodes(:,2); % Inlet X velocity
% bc_f.drichletVector(3*mesh_f.inletNodes-1) = fluidInletNodes(:,3); % Inlet Y velocity
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
% mesh_f.fsiDOF = [2*mesh_f.fsiNodes-1; 2*mesh_f.fsiNodes];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% ALE MESH  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basically all the nodes and elements are same as the fluid mesh.
mesh_ale    = mesh_f;
bc_ale      = bc_f;

%% Formulating the FSI nodes DOF
bc_ale.fsiDOF = zeros(2*length(mesh_ale.fsiNodes),1);
bc_ale.fsiDOF(1:2:end) = 2*mesh_ale.fsiNodes-1;
bc_ale.fsiDOF(2:2:end) = 2*mesh_ale.fsiNodes;
nodesOnCyl = [1793,1801,1813,1828,1847,1859,1877,1895,1912,1926,1946,1959,1978,1987,1990,1984,1974,1956,1940,1922,1905,1892,1874,1857,1846,1829,1815,1804];
for i=1:length(nodesOnCyl)
    tempNode = find(mesh_ale.globalWallNumbers == nodesOnCyl(i));
    localNodesOnCyl(i) = mesh_ale.wallNodes(tempNode);
    mesh_ale.wallNodes(tempNode) = [];
    mesh_ale.globalWallNumbers(tempNode) = [];
end
localNodesOnCyl = localNodesOnCyl';
% bc_ale.fsiDOF = [2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes];

%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
% Here we constarin the X movement on the inlet and outlet, Y movement on
% the top and bottom walls. This way we prevent bad mesh movement.
% bc_ale.drichletDOF = [2*mesh_ale.inletNodes-1; 2*mesh_ale.inletNodes; 2*mesh_ale.wallNodes-1; 2*mesh_ale.wallNodes; 2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes; 2*mesh_ale.outletNodes-1; 2*mesh_ale.outletNodes];
bc_ale.drichletDOF = [2*mesh_ale.inletNodes-1; 2*mesh_ale.wallNodes; 2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes; 2*mesh_ale.outletNodes-1;2*localNodesOnCyl;2*localNodesOnCyl-1];

% Drichlet values vector
bc_ale.drichletVector = zeros(2*length(mesh_ale.nodes),1);

% Finding the DOFs which are to be computed
bc_ale.freeDOF = setdiff(1:2*length(mesh_ale.nodes), bc_ale.drichletDOF);

%% ::: IMPORTANT :: IMPORTANT :: 
% The FSI interface nodes Drichlet vector should always be updated in every
% time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %% Formulating entities for Structure mesh
% % % mesh_s.nodes = structureNodes(:,2:4);
% % % mesh_s.nodesGlobalNumbers  = structureNodes(:,1);
% % % mesh_s.elements = structureElements(:,2:4);
% % % mesh_s.fixedNodes = structureFixedNodes(:,1);
% % % mesh_s.forceNodes = structureForceNodes(:,1);
% % % mesh_s.fsiNodes = fsiInterfaceNodes(:,1);
% % % mesh_s.fsiNodesGlobal = mesh_s.fsiNodes;


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
% bc_s.neumannVector(2*mesh_s.forceNodes-1) = structureForceNodes(:,2); % X Component of force
% bc_s.neumannVector(2*mesh_s.forceNodes) = structureForceNodes(:,3);   % Y Component of force

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

fprintf('\nTotal Fluid Nodes :: %d \n',length(mesh_f.nodes));
fprintf('Total Fluid Elements :: %d \n', length(mesh_f.elements));

fprintf('\nTotal Stricture Nodes :: %d \n',length(mesh_s.nodes));
fprintf('Total Structure Elements :: %d \n', length(mesh_s.elements));


% End of the method
end