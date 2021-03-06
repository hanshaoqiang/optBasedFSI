function [mesh_s, bc_s, physics_s, mesh_f, physics_f, bc_f, mesh_ale, bc_ale] = readMeshFromGIDOutput(casename)
%% This function reads the output files form the GID output and formulates
% the necessary arrays and element tables for both fluid and structural
% solvers.
% 
% The mesh files are by default assumed to be in the path 
% /mesh/GID_output
%   INPUT 
%       casename    : name of the gid case name where mesh is generated.
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

meshFileName = [casename '.dat'];
fstring = fileread(meshFileName); % read the mesh file as one string

%% FLUID_DENSITY
fblocks = regexp(fstring,'FLUID_DENSITY','split');
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.rho = out(:,1);

%% FLUID_DYNAMIC_VISCOSITY
fblocks = regexp(fstring,'FLUID_DYNAMIC_VISCOSITY','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.mu = out(:,1);
physics_f.nu = physics_f.mu/physics_f.rho; % Kinematic Viscosity

%% FLUID_SOLVER
fblocks = regexp(fstring,'FLUID_SOLVER','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.transient_f.solver = out(:,1);

%% FLUID_TIME_INTEGRATION
fblocks = regexp(fstring,'FLUID_TIME_INTEGRATION','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.transient_f.integrationScheme = out(:,1);

%% FLUID_START_TIME
fblocks = regexp(fstring,'FLUID_START_TIME','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.transient_f.Tstart = out(:,1);


%% FLUID_END_TIME
fblocks = regexp(fstring,'FLUID_END_TIME','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.transient_f.Tend = out(:,1);


%% FLUID_TIME_STEP
fblocks = regexp(fstring,'FLUID_TIME_STEP','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_f.transient_f.dt = out(:,1);


%% FLUID_NODES
fblocks = regexp(fstring,'FLUID_NODES','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
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
mesh_f.elements = out(:,2:4);


%% FLUID_DIRICHLET_CONDITIONS
fblocks = regexp(fstring,'FLUID_DIRICHLET_CONDITIONS','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_f.drichletNodes = out(:,1);
mesh_f.drichletConditions = out(:,2:4);

%% FLUID_NEUMANN_CONDITIONS
fblocks = regexp(fstring,'FLUID_NEUMANN_CONDITIONS','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_f.neumannNodes = out(:,1);
mesh_f.neumannConditions = out(:,2:4);


%% STRUCTURE_DENSITY
fblocks = regexp(fstring,'STRUCTURE_DENSITY','split');
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.rho = out(:,1);

%% YOUNGS_MODULUS
fblocks = regexp(fstring,'YOUNGS_MODULUS','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.E = out(:,1);

%% POISSON_RATIO
fblocks = regexp(fstring,'POISSON_RATIO','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.nu = out(:,1);

%% STRUCTURAL_SOLVER
fblocks = regexp(fstring,'STRUCTURAL_SOLVER','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.transient_s.solver = out(:,1);

%% STRUCTURAL_TIME_INTEGRATION
fblocks = regexp(fstring,'STRUCTURAL_TIME_INTEGRATION','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.transient_s.integrationScheme = out(:,1);

%% STRUCTURAL_START_TIME
fblocks = regexp(fstring,'STRUCTURAL_START_TIME','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.transient_s.Tstart = out(:,1);


%% STRUCTURAL_END_TIME
fblocks = regexp(fstring,'STRUCTURAL_END_TIME','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.transient_s.Tend = out(:,1);


%% STRUCTURAL_TIME_STEP
fblocks = regexp(fstring,'STRUCTURAL_TIME_STEP','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',',','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
physics_s.transient_s.dt = out(:,1);


%% STRUCTURAL_NODES
fblocks = regexp(fstring,'STRUCTURAL_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_s.nodes = out(:,2:4);
mesh_s.nodesGlobalNumbers = out(:,1);


%% STRUCTURAL_ELEMENTS
fblocks = regexp(fstring,'STRUCTURAL_ELEMENTS','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_s.elements = out(:,2:4);


%% STRUCTURAL_DIRICHLET_CONDITIONS
fblocks = regexp(fstring,'STRUCTURAL_DIRICHLET_CONDITIONS','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_s.drichletNodes = out(:,1);
mesh_s.drichletConditions = out(:,2:4);

%% STRUCTURAL_NEUMANN_CONDITIONS
fblocks = regexp(fstring,'STRUCTURAL_NEUMANN_CONDITIONS','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_s.neumannNodes = out(:,1);
mesh_s.neumannConditions = out(:,2:4);

%% FSI_INTERFACE_NODES
fblocks = regexp(fstring,'FSI_INTERFACE_NODES','split'); 
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
mesh_s.fsiNodesGlobal = out(:,1);
mesh_f.fsiNodesGlobal = out(:,1);
mesh_s.fsiNodes = out(:,1);
mesh_f.fsiNodes = out(:,1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% FLUID %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% For Drichlet Nodes
for i = 1:length(mesh_f.drichletNodes)
    node_global = mesh_f.drichletNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.drichletNodes(i) = node_local;
end

% For Neumann Nodes
for i = 1:length(mesh_f.neumannNodes)
    node_global = mesh_f.neumannNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.neumannNodes(i) = node_local;
end

% For fsiNodes
for i = 1:length(mesh_f.fsiNodes)
    node_global = mesh_f.fsiNodes(i);
    node_local = find(node_global == mesh_f.nodesGlobalNumbers);
    mesh_f.fsiNodes(i) = node_local;
end


%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
bc_f.drichletDOF = [3*mesh_f.drichletNodes-2; 3*mesh_f.drichletNodes-1];

% Drichlet values vector
bc_f.drichletVector = zeros(3*length(mesh_f.nodes),1);

maxVel = 2*norm(fluidInletNodes(:,2),inf);
y = mesh_f.nodes(mesh_f.inletNodes,2);
parabolicVelProf = 1.5*maxVel*(4.0/0.1681)*y.*(0.41-y);

bc_f.drichletVector(3*mesh_f.inletNodes-2) = mesh_f.drichletConditions(:,1); % X Condition
bc_f.drichletVector(3*mesh_f.inletNodes-1) = mesh_f.drichletConditions(:,2); % Y Condition

% IMPORTANT :: IMPORTANT :: Keep in mind that the fsi Values should be
% updated at every time step


% Neumann DOF vector
bc_f.neumannDOF = 3*mesh_f.neumannNodes;

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
% bc_ale.fsiDOF = [2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes];

%% Formulating the drichletDOFs and neumannDOFs vector and corresponding full values vectors
% Drichlet DOF vector
% Here we constarin the X movement on the inlet and outlet, Y movement on
% the top and bottom walls. This way we prevent bad mesh movement.
% bc_ale.drichletDOF = [2*mesh_ale.inletNodes-1; 2*mesh_ale.inletNodes; 2*mesh_ale.wallNodes-1; 2*mesh_ale.wallNodes; 2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes; 2*mesh_ale.outletNodes-1; 2*mesh_ale.outletNodes];
bc_ale.drichletDOF = [2*mesh_ale.inletNodes-1; 2*mesh_ale.wallNodes; 2*mesh_ale.fsiNodes-1; 2*mesh_ale.fsiNodes; 2*mesh_ale.outletNodes-1];

% Drichlet values vector
bc_ale.drichletVector = zeros(2*length(mesh_ale.nodes),1);

% Finding the DOFs which are to be computed
bc_ale.freeDOF = setdiff(1:2*length(mesh_ale.nodes), bc_ale.drichletDOF);

%% ::: IMPORTANT :: IMPORTANT :: 
% The FSI interface nodes Drichlet vector should always be updated in every
% time step












end