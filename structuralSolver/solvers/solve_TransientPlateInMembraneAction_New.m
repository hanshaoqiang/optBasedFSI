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
function [displacement, velocity, acceleration] = solve_TransientPlateInMembraneAction_New(mesh, bc, physics, transient, dispPrevious, velPrevious, accPrevious, fFSI, fFSIprev)
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
F = bc.neumannVector; % Force Vector 
Fprev = bc.neumannVector; % Force Vector 

%% Applying the FSI Interface force on the Force vector
F(2*mesh.fsiNodes-1) = F( 2*mesh.fsiNodes-1) + fFSI(1:2:end);
F(2*mesh.fsiNodes) = F( 2*mesh.fsiNodes) + fFSI(2:2:end);

%% Applying the FSI Interface force on the Force vector
Fprev(2*mesh.fsiNodes-1) = Fprev(2*mesh.fsiNodes-1) + fFSIprev(1:2:end);
Fprev(2*mesh.fsiNodes) = Fprev(2*mesh.fsiNodes) + fFSIprev(2:2:end);


%% 1. Compute the master stiffness matrix of the structure
[K, M] = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,physics,physics);


% System Matrix
G = (6/dt/dt)*M + K;
% Right hand side vector
RHS = F + Fprev -K*dispPrevious + M*accPrevious + (6/dt/dt)*M*(dispPrevious + dt*velPrevious);

%% Applying the Drichlet boundary conditions on acceleration
RHS = RHS - G * bc.drichletVector;

%% 3. Compute acceleration vector
displacement(bc.freeDOF) = G(bc.freeDOF,bc.freeDOF)\RHS(bc.freeDOF);


%% 4. Assemble to the complete displacement vector
displacement(bc.drichletDOF) = bc.drichletVector(bc.drichletDOF);
displacement = displacement';

%% Computing the acceleration
acceleration = (6/dt/dt)*(displacement -dispPrevious - dt*velPrevious - (dt*dt/3)*accPrevious);


%% Computing the velocity vector
velocity = velPrevious + dt*0.5*(accPrevious + acceleration);


end
