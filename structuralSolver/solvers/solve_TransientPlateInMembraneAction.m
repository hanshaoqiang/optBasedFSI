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
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement, velocity, acceleration,sens_fsi,Vel_new] = solve_TransientPlateInMembraneAction(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious, fFSI)
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
%        accPrevious : acceleration in previous time step
%               fFSI : force vector on the FSI interface a calculated from
%                      the fluid solver
%           fFSIprev : force vector on the FSI interface in previous time step a calculated from
%                      the fluid solver
%
%             OUTPUT
%       displacement : The resulting displacement field at current time
%           velocity : The resulting velocity field at current time
%
%% Function main body

%% Setting up initial values required for simulation
t = transient.Tstart;
dt = transient.dt;
F = [zeros(2*length(mesh.nodes),1); bc.neumannVector]; % Force Vector 

%% Applying the FSI Interface force on the Force vector
F(2*length(mesh.nodes) + (2*mesh.fsiNodes-1)) = F(2*length(mesh.nodes) + (2*mesh.fsiNodes-1)) + fFSI(1:2:end);
F(2*length(mesh.nodes) + 2*mesh.fsiNodes) = F(2*length(mesh.nodes) + 2*mesh.fsiNodes) + fFSI(2:2:end);


%% Formulating previous U
Uprev = [dispPrevious; velPrevious];

%% Formulating previous Udot
UdotPrev = [velPrevious; accPrevious];


%% 1. Compute the master stiffness matrix of the structure
[K, M, bfVector] = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,physics,physics);
M = 0.5*M;
% Adding the body force contribution to the force vector
F(2*length(mesh.nodes)+1:end) = F(2*length(mesh.nodes)+1:end) + bfVector; 

%% Formulating the matrices(simultanious equations) for time integration
P = [eye(2*length(mesh.nodes)) zeros(2*length(mesh.nodes));
     zeros(2*length(mesh.nodes)) M];
 
Q = [zeros(2*length(mesh.nodes)) -eye(2*length(mesh.nodes));
     K zeros(2*length(mesh.nodes))];

%% Reformulaing the boundary conditions to include the conditions for momentum term also
bc.drichletVector = [bc.drichletVector; zeros(2*length(mesh.nodes),1) ];
bc.drichletDOF = [bc.drichletDOF; 2*length(mesh.nodes) + bc.drichletDOF];
bc.freeDOF = [bc.freeDOF, 2*length(mesh.nodes) + bc.freeDOF];

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

%% Calculating the senstitivities dUdot_dU
% % % % % % % % Extracting the velocity part of G
% % % % % % %  dispDOF  = 1:2*length(mesh.nodes);
% % % % % % %  velDOF = 2*length(mesh.nodes)+1 : 4*length(mesh.nodes);
% % % % % % % % G_vel_vel   = G(velDOF,velDOF);
% % % % % % % % G_vel_disp  = G(velDOF,dispDOF);
% % % % % % % G_disp_vel  = G(dispDOF,velDOF);
% % % % % % % G_disp_disp = G(dispDOF,dispDOF);
% % % % % % % % dUdot_dU_2  = -1*(G_vel_vel^-1)*G_vel_disp;
% % % % % % % dUdot_dU    = -1*(G_disp_vel^-1)*G_disp_disp;
% % % % % % % dUdot_dU_fsi    = dUdot_dU(mesh.fsiDOF, mesh.fsiDOF);
% % % % % % % % % % % DeMat = zeros(4*length(mesh.nodes),2*length(mesh.nodes));
% % % % % % % % % % % x = zeros(4*length(mesh.nodes),2*length(mesh.nodes));
% % % % % % % % % % % % Calculating the RHS 
% % % % % % % % % % % DeMat(1:2*length(mesh.nodes),:) = eye(2*length(mesh.nodes));
% % % % % % % % % % % RHS_ses = -G*DeMat;
% % % % % % % % % % % % Solving for sensitivities
% % % % % % % % % % % x(2*length(mesh.nodes)+1:4*length(mesh.nodes),:) = G(2*length(mesh.nodes)+1:4*length(mesh.nodes),2*length(mesh.nodes)+1:4*length(mesh.nodes))\RHS_ses(2*length(mesh.nodes)+1:4*length(mesh.nodes),:);
dispDOF  = 1:2*length(mesh.nodes);
velDOF = 2*length(mesh.nodes)+1 : 4*length(mesh.nodes);
G_vel_vel   = G(velDOF,velDOF);
G_vel_disp  = G(velDOF,dispDOF);
G_disp_vel  = G(dispDOF,velDOF);
G_disp_disp = G(dispDOF,dispDOF);

L = [G_disp_vel zeros(2*length(mesh.nodes));
     G_vel_vel  -eye(2*length(mesh.nodes))];
 
RHS_ses = -1*[G_disp_disp*eye(2*length(mesh.nodes));
              G_vel_disp*eye(2*length(mesh.nodes))];
          
% Solving for sensitivities
sens = L\RHS_ses;
sens = sens(1:2*length(mesh.nodes),:);
sens_fsi = sens(mesh.fsiDOF,mesh.fsiDOF);

%%% Verifying the above with finite difference.
delta = 1E-4;
% Petrubing Disp on Drichlet BC
nodeDOF = 2*mesh.fsiNodes(10);
% Solving for new velocity.
% Extracting the displacement vector from Drichlet Vector
U_new = U(1:2*length(mesh.nodes));
Vel_old = U(2*length(mesh.nodes)+1:4*length(mesh.nodes));
U_new(nodeDOF) = U_new(nodeDOF) + delta;

% Using the changed displacement to calculate the new velocity
Vel_new = -G_disp_vel\(G_disp_disp*U_new);

% Calculating diff
diff = Vel_new - Vel_old;

dUdot_dU_fd = diff./delta;
dUdot_dU_3 = full(sens(:,nodeDOF));

% Extracting only the fsi velocity
Vel_new = Vel_new(mesh.fsiDOF);
end
