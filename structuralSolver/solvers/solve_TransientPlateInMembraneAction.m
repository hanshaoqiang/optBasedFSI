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
function [displacement, velocity, acceleration] = solve_TransientPlateInMembraneAction(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious, fFSI)
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
displacement = [];
msgEmpty = '';
reverseStr = '';
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

end
