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
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de)         %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement, velocity, acceleration] = solve_MeshMovement(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious)
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
U = zeros(4*length(mesh.nodes),1);
t = transient.Tstart;
dt = transient.dt;
F = zeros(4*length(mesh.nodes),1); % Force Vector 

%% Formulating the Drichlet boundary condition vector
bc.drichletVector = [bc.drichletVector; zeros(2*length(mesh.nodes),1) ];
bc.drichletVector(bc.fsiDOF) = bc.fsiDisp;
bc.drichletVector(2*length(mesh.nodes) + bc.fsiDOF) = bc.fsiVel;
% bc.drichletDOF = [bc.drichletDOF; 2*length(mesh.nodes) + bc.drichletDOF];
% bc.freeDOF  = [bc.freeDOF, 2*length(mesh.nodes) + bc.freeDOF];
bc.freeDOF  = setdiff(1:4*length(mesh.nodes), bc.drichletDOF);

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


%% Computing the velocity
Udot = a*(U - Uprev) - d*UdotPrev;

%% 5. compute the complete displacement, velocity, acceleration vectors
displacement = U(1:2*length(mesh.nodes));
velocity = U(2*length(mesh.nodes)+1:end);
velocityTest = Udot(1:2*length(mesh.nodes));
acceleration = Udot(2*length(mesh.nodes)+1:end);

end
