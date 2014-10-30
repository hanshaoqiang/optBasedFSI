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
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de) 
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement] = solve_steadyStateStructure(mesh,bc,F,materialProperties,analysis)
%% Function documentation
%
% Returns the displacement field corresponding to a plain stress/strain
% analysis for the given mesh of the geometry together with its Dirichlet
% and Neumann boundary conditions.
%
%              Input :
%               mesh : Elements and nodes of the mesh
%                 bc : Boundary conditions formulatied in terms of nodal
%                       DOFs
%                  F : Global load vector(obtained from the fluid)
% materialProperties : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%      
%             Output :
%       displacement : The resulting displacement field
%
%
%% Function main body

%% 1. Compute the master stiffness matrix of the structure
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

% Assigning the RHS of the system
RHS = F;

%% Applying the Drichlet boundary conditions
RHS = RHS - K * bc.drichletVector;

% Solving the system
displacement(bc.freeDOF) = K(bc.freeDOF,bc.freeDOF)\RHS(bc.freeDOF);

% Assembling the solution vector
displacement(bc.drichletDOF) = bc.drichletVector(bc.drichletDOF);

%% 5. compute the complete load vector and verify the results
% Compute the complete force vector
Fcomplete = K*displacement';

% Verify the resutls
F_verification = Fcomplete;

%  A tolerance value
tolerance = 1e6;

% Compute a residual
residual = norm(F_verification-F);

if residual>=tolerance
   fprintf('\t The computed load vector has been detected to be different \n');
   fprintf('\t than the global load vector in the range of %d \n',residual);
end

end