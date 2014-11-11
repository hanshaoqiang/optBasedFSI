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
elseDOF = setdiff(1:2*length(mesh.nodes),fsiDOF);


% % % % fsiDOF_disp = bc.fsiDOF;
% % % % 
% % % % elseDOF_disp= setdiff(1:2*length(mesh.nodes), fsiDOF_disp);
% % % % fsiDOF_vel  = 2*length(mesh.nodes) + bc.fsiDOF;
% % % % elseDOF_vel = setdiff(2*length(mesh.nodes)+1:4*length(mesh.nodes), fsiDOF_vel);

% % % % % % % Extracting corresponding parts from G matrix
% % % % % % G_dede      = full(G(elseDOF_disp,elseDOF_disp));
% % % % % % G_didi      = full(G(fsiDOF_disp,fsiDOF_disp));
% % % % % % G_dpedpe    = full(G(elseDOF_vel, elseDOF_vel));
% % % % % % G_dpidpi    = full(G(fsiDOF_vel,fsiDOF_vel));
% % % % % % G_dedi      = full(G(elseDOF_disp,fsiDOF_disp));
% % % % % % G_dide      = full(G(fsiDOF_disp,elseDOF_disp));
% % % % % % G_dpedpi    = full(G(elseDOF_vel,fsiDOF_vel));
% % % % % % G_dpidpe    = full(G(fsiDOF_vel,elseDOF_vel));
% % % % % % G_dpedi     = full(G(elseDOF_vel,fsiDOF_disp));
% % % % % % G_dpidi     = full(G(fsiDOF_vel,fsiDOF_disp));
% % % % % % 
% % % % % % L = [G_dpedpe G_dpedpi;
% % % % % %      G_dpidpe G_dpidpi ];
% % % % % % I = eye(length(bc.fsiDOF));
% % % % % % 
% % % % % % R1 = -1*[G_dpedi*I;         G_dpidi*I];
% % % % % % R2 = -1*[G_dpedpi*duip_dui; G_dpidpi*duip_dui];
% % % % % % 
%B = [G_dede G_dedi; G_dide G_didi];% % % % % % R = R1+R2;
% % % % % % mesh_ses = (L^-1) * R;
%% Full analysis with disp and velocity
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses = zeros(4*length(mesh.nodes),length(fsiDOF_disp));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % DeMat = zeros(4*length(mesh.nodes), length(fsiDOF_disp));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % DeMat(fsiDOF_disp,:) = eye(length(fsiDOF_disp));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % DeMat(fsiDOF_vel,:)  = duip_dui;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % RHS_ses = -G*DeMat;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Solving for sensitivity
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % free_ses_dof = [elseDOF_disp,elseDOF_vel];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % x = G(free_ses_dof,free_ses_dof)\RHS_ses(free_ses_dof,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses(free_ses_dof,:) = x;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Putting bc again in
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses(fsiDOF_disp,:) = eye(length(fsiDOF_disp));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses(fsiDOF_vel,:)  = duip_dui;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses_1 = mesh_ses(elseDOF_disp,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % fsiNode = 10;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh_ses_3  = full(mesh_ses_1(:,2*fsiNode));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % delta = 10^(-4);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % bc_new = bc;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % U_new = zeros(4*length(mesh.nodes),1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % RHS_new = F - P*(b*UdotPrev - c*Uprev);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %FSI DOF to pertrub.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % pet_node = mesh.fsiNodes(fsiNode);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % pet_DOF = 2*pet_node;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %Petrubing
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % bc_new.drichletVector(pet_DOF) = bc_new.drichletVector(pet_DOF) + delta;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Apply the velocity with pretrubed displacements as Drichlet velocity.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % bc_new.drichletVector(2*length(mesh.nodes) + bc.fsiDOF) = petVel;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %Applying the Drichlet boundary conditions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % RHS_new = RHS_new - G * bc_new.drichletVector;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %3. Compute displacement vector
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % U_new(bc_new.freeDOF) = G(bc_new.freeDOF,bc_new.freeDOF)\RHS_new(bc_new.freeDOF);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %4. Assemble to the complete displacement vector
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % U_new(bc_new.drichletDOF) = bc_new.drichletVector(bc_new.drichletDOF);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %Calcualting the gradient
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % y = (U_new - U)./delta;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % y = y(elseDOF_disp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % err = norm(y - mesh_ses_3, inf)
%% Analysis by neglecting velocity
% % G_dede      = K(elseDOF,elseDOF);
% % G_dedi      = K(elseDOF,fsiDOF);
% % G_didi      = K(fsiDOF, fsiDOF);
% % G_dide      = K(fsiDOF, fsiDOF);
% % mesh_ses    = -1*G_dede\G_dedi;
UK = zeros(2*length(mesh.nodes),1);
UK(2*mesh.fsiNodes(20),1) = 1;
rhs = -K*UK;
mesh_ses= zeros(2*length(mesh.nodes),1);
mesh_ses(elseDOF) = K(elseDOF,elseDOF)\rhs(elseDOF);
mesh_ses(2*mesh.fsiNodes(20),:) = 1;


% Finite difference
delta = 1E-5;
pet_node = mesh.fsiNodes(20);
pet_DOF = 2*pet_node;
B = K;
bc_bu.drichletVector(bc_bu.fsiDOF) = bc_bu.fsiDisp;
dv = zeros(2*length(mesh.nodes),1);
dv(2*mesh.fsiNodes(20),1) = 10000000;
rhs = -B*dv;
% Solving
displac(bc_bu.freeDOF) = B(bc_bu.freeDOF,bc_bu.freeDOF)\rhs(bc_bu.freeDOF);
displac(1,2*mesh.fsiNodes(20)) = 10000000;

% Perturbing the fsiNode
dv(2*mesh.fsiNodes(20),1) = dv(2*mesh.fsiNodes(20),1) + delta;
bc_bu.drichletVector(pet_DOF) = bc_bu.drichletVector(pet_DOF) + delta;
rhs = -B*dv;
% Solving
displac_pert(bc_bu.freeDOF) = B(bc_bu.freeDOF,bc_bu.freeDOF)\rhs(bc_bu.freeDOF);
displac_pert(1,2*mesh.fsiNodes(20)) = dv(2*mesh.fsiNodes(20),1) + delta;

mesh_sens_fd = (displac_pert - displac)./delta;
mesh_sens_fd = mesh_sens_fd';
mesh_sens_ana = mesh_ses(:,:);


end
