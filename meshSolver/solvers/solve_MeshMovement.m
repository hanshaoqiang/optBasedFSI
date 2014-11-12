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
%   Aditya Ghantasala M.Sc             (aditya.ghantasala@tum.de)         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement, velocity, acceleration] = solve_MeshMovement(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious,duip_dui,petVel)
%% Function documentation
%
% Returns the displacement field corresponding to a plain stress/strain
% analysis for the given mesh of the geometry together with its Dirichlet
% and Neumann boundary conditions.
%
%              INPUT 
%               mesh : Elements and nodes of the mesh
%                 bc : structure containing the boundary conditions applied
%                      interms of node numbering
%            physics : The material properties of the structure
%          transient : The properties time integration
%       dispPrevious : displacement from previous time step
%        velPrevious : velocity from previous time step
%        accPrevious : Acceleration from previous time step
%           duip_dui : Sensitivity of the structural solver over the FSI
%                      interface.
%               
%
%             OUTPUT
%       displacement : The resulting displacement field at current time
%           velocity : The resulting velocity field at current time
%       acceleration : The resulting acceleration field at current time.
%
%% Function main body

%% Setting up initial values required for simulation
displacement = [];
t = transient.Tstart;
dt = transient.dt;
F = zeros(4*length(mesh.nodes),1); % Force Vector 
bc_bu = bc;

%% Formulating the Drichlet boundary condition vector
bc.drichletVector = [bc.drichletVector; zeros(2*length(mesh.nodes),1) ];
bc.drichletVector(bc.fsiDOF) = bc.fsiDisp;
bc.drichletVector(2*length(mesh.nodes) + bc.fsiDOF) = bc.fsiVel;
bc.drichletDOF = [bc.drichletDOF; 2*length(mesh.nodes) + bc.drichletDOF];
bc.freeDOF = [bc.freeDOF, 2*length(mesh.nodes) + bc.freeDOF];

%% Formulating previous U
Uprev = [dispPrevious; velPrevious];

%% Formulating previous Udot
UdotPrev = [velPrevious; accPrevious];


%% 1. Compute the master stiffness matrix of the structure
[K, M] = computeStiffnessMatrixMeshMotionLinear(mesh,physics,physics);


%% Formulating the matrices(simultanious equations) for time integration
P = [eye(2*length(mesh.nodes)) zeros(2*length(mesh.nodes));
     zeros(2*length(mesh.nodes)) M];
 
Q = [zeros(2*length(mesh.nodes)) -eye(2*length(mesh.nodes));
     K zeros(2*length(mesh.nodes))];

%% Calculating the Bossak time integration specific matrices
a = 1/(transient.gamma*transient.dt);
b = transient.alphaBeta  - ((1-transient.gamma)*(1-transient.alphaBeta))/transient.gamma;
c = (1-transient.alphaBeta)*a;
d = (1-transient.gamma)/transient.gamma;

% System Matrix
G = c*P + Q;
% Right hand side vector
RHS = F - P*(b*UdotPrev - c*Uprev);

%% Applying the Drichlet boundary conditions
RHS = RHS - G * bc.drichletVector;

%% 3. Compute displacement vector
U(bc.freeDOF) = G(bc.freeDOF,bc.freeDOF)\RHS(bc.freeDOF);


%% 4. Assemble to the complete displacement vector
U(bc.drichletDOF) = bc.drichletVector(bc.drichletDOF);
U = U';

%% Computing the velocity
Udot = a*(U - Uprev) - d*UdotPrev;

%% 5. compute the complete displacement, velocity, acceleration vectors
displacement = U(1:2*length(mesh.nodes));
velocity = U(2*length(mesh.nodes)+1:end);
velocityTest = Udot(1:2*length(mesh.nodes));
acceleration = Udot(2*length(mesh.nodes)+1:end);


%% Writing Fluid mesh displacement Result
mesh.elemOrder = 3;
FluidMeshDispTitle = ['Fluid_Mesh_Disp_Result_at_Time_',num2str(transient.i),'_',num2str(0),'.vtk'];
formulatedFMeshResult = zeros(3*length(mesh.nodes),1);
formulatedFMeshResult(3*(1:length(mesh.nodes))-2) = displacement(1:2:end);
formulatedFMeshResult(3*(1:length(mesh.nodes))-1) = displacement(2:2:end);
formulatedFMeshResult(3*mesh.fsiNodes-2) = 1;
formulatedFMeshResult(3*mesh.fsiNodes-1) = 1;

writeResultToVtk(mesh, formulatedFMeshResult, FluidMeshDispTitle);
clear formulatedFMeshResult



%% Formulating for term dd/dd_fsi
fsiDOF = zeros(2*length(mesh.fsiNodes),1);
fsiDOF(1:2:end) = 2*mesh.fsiNodes-1;
fsiDOF(2:2:end) = 2*mesh.fsiNodes;
elseDOF = setdiff(1:2*length(mesh.nodes),fsiDOF); % This is not correct because it should not contain the fixed wall DOFs

%% Analysis by neglecting velocity
UK = zeros(2*length(mesh.nodes),2*length(mesh.fsiNodes));
UK(fsiDOF,:) = eye(length(fsiDOF));
rhs = -K*UK;
mesh_ses = zeros(2*length(mesh.nodes),2*length(mesh.fsiNodes));
mesh_ses(bc_bu.freeDOF,:) = K(bc_bu.freeDOF,bc_bu.freeDOF)\rhs(bc_bu.freeDOF,:);
mesh_ses(fsiDOF,:) = eye(length(fsiDOF));


% Finite difference
delta = 1E-4;
pet_node = mesh.fsiNodes(20);
pet_DOF = 2*pet_node;
B = K;
bc_bu.drichletVector(bc_bu.fsiDOF) = bc_bu.fsiDisp;
rhs = -B*bc_bu.drichletVector;
% Solving
displac(bc_bu.freeDOF) = B(bc_bu.freeDOF,bc_bu.freeDOF)\rhs(bc_bu.freeDOF);
displac(bc_bu.drichletDOF) = bc_bu.drichletVector(bc_bu.drichletDOF);

% Perturbing the fsiNode
bc_bu.drichletVector(pet_DOF) = bc_bu.drichletVector(pet_DOF) + delta;
rhs = -B*bc_bu.drichletVector;
% Solving
displac_pert(bc_bu.freeDOF) = B(bc_bu.freeDOF,bc_bu.freeDOF)\rhs(bc_bu.freeDOF);
displac_pert(bc_bu.drichletDOF) = bc_bu.drichletVector(bc_bu.drichletDOF);

mesh_sens_fd = (displac_pert - displac)./delta;
mesh_sens_fd = mesh_sens_fd';
mesh_sens_ana = mesh_ses(:,40);


end
