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
function [displacement, velocity, acceleration] = solve_MeshMovement(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious,duip_dui)
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
%               
%
%             OUTPUT
%       displacement : The resulting displacement field at current time
%           velocity : The resulting velocity field at current time
%
%% Function main body

%% Setting up initial values required for simulation
displacement = [];
t = transient.Tstart;
dt = transient.dt;
F = zeros(4*length(mesh.nodes),1); % Force Vector 

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


%% Formulating for term dd/dd_fsi
fsiDOF_disp = bc.fsiDOF;
elseDOF_disp= setdiff(1:2*length(mesh.nodes), fsiDOF_disp);
fsiDOF_vel  = 2*length(mesh.nodes) + bc.fsiDOF;
elseDOF_vel = setdiff(2*length(mesh.nodes)+1:4*length(mesh.nodes), fsiDOF_vel);

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
% % % % % % R = R1+R2;
% % % % % % mesh_ses = (L^-1) * R;
mesh_ses = zeros(4*length(mesh.nodes),length(fsiDOF_disp));
DeMat = zeros(4*length(mesh.nodes), length(fsiDOF_disp));
DeMat(fsiDOF_disp,:) = eye(length(fsiDOF_disp));
DeMat(fsiDOF_vel,:)  = duip_dui;
RHS_ses = -G*DeMat;

% Solving for sensitivity
free_ses_dof = [elseDOF_disp,elseDOF_vel];
x = G(free_ses_dof,free_ses_dof)\RHS_ses(free_ses_dof,:);
mesh_ses(free_ses_dof,:) = x;
% Putting bc again in
% mesh_ses(fsiDOF_disp,:) = eye(length(fsiDOF_disp));
% mesh_ses(fsiDOF_vel,:)  = duip_dui;
mesh_ses_1 = mesh_ses(elseDOF_disp,:);


%%% Checking the above formulation with finite difference.
fsiNode = 10;
mesh_ses_3  = full(mesh_ses_1(:,2*fsiNode));
err = 0;
for p = 1 : 20
    
    delta = 10^(-p);
    bc_new = bc;
    U_new = zeros(4*length(mesh.nodes),1);
    RHS_new = F - P*(b*UdotPrev - c*Uprev);
    % FSI DOF to pertrub.
    pet_node = mesh.fsiNodes(fsiNode);
    pet_DOF = 2*pet_node;
    % Petrubing
    bc_new.drichletVector(pet_DOF) = bc_new.drichletVector(pet_DOF) + delta;
    
    % Applying the Drichlet boundary conditions
    RHS_new = RHS_new - G * bc_new.drichletVector;
    
    % 3. Compute displacement vector
    U_new(bc_new.freeDOF) = G(bc_new.freeDOF,bc_new.freeDOF)\RHS_new(bc_new.freeDOF);
    
    
    % 4. Assemble to the complete displacement vector
    U_new(bc_new.drichletDOF) = bc_new.drichletVector(bc_new.drichletDOF);
    %U_new = U_new';
    % Calcualting the gradient
    y = (U_new - U)./delta;
    mesh_ses_fd(:,p) = y(elseDOF_disp);
    
    err(p) = norm(mesh_ses_fd(:,p) - mesh_ses_3, inf);
end







end
