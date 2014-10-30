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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [up, fFSI] = solve_SteadyStateNavierStokesProblem2D(mesh, physics, bc, upInitial)
%% Function documentation
%
% Solve the steady state Navier-Stokes problem in 2D using triangular basis
% basis for the spatial discretization.
%
%   INPUT:
%       geometry                :   Structure containing the geometry
%                                   specified through NURBS
%       mesh                    :   Strucure containing the mesh(nodes, elements, boundaryElements ...)
%       physics                 :   Sturcure containing the physicsl
%                                   properties of the simulation setup. (also includes the mesh size)
%       bc                      :   Strucure conaining the specified
%                                   Boundary conditions.
%
%   OUTPUT:
%       up                      :   Solution of the Steady state problem
%       fFSI                    :   Force on the FSI interface.
%   
%
%% Function main body

%% Read input
nCPs = length(mesh.nodes);

% Compute the number of degrees of freedom
nDoFs = 3*nCPs;

% Initialize the global flux vector
F = zeros(nDoFs,1);


%% Get the initial values for the discrete solution vector and its rate

% Get the initial condition for the scalar concetration field
up = upInitial;

%% Preamble

% Initialize string for the Newton-Rapson iterations
reverseStrNR = '';


%% Loop over all the Newton-Rapson iterations
msgPNR = sprintf('\t \t Looping over all the Newton iterations \n \t \t -------------------------------------- \n \n');
fprintf(msgPNR);
for i_newton = 1:physics.nonLinearScheme.maxIter
    %% Compute the tangent stiffness, the stiffness and the mass matrix of the finite element system
    % Compute the stiffness, tangent stiffness, mass matrix and body
    % force vector at the current nonlinear iteration step
    [K,KTangent,F] = computeSteadyStateMatricesForIncompressibleNavierStokes2D(up,mesh, physics);
    
    % Compute the left-hand side nonlinear stiffness matrix for the
    % correction step
    G = K+KTangent;
    
    %% Updating the right hand side vector to account for Drichlet and Neumann boundary conditons
    r = K*up - F;
    
    RHS = - r;
    
    % Compute the generalized right-hand side load vector
    if (i_newton==1)
        % Compute the right-hand side load vector due to the inhomogeneous boundary
        % conditions
        % Update the right-hand side vector using the prescribed values
        RHS = RHS - G*bc.drichletVector;
    end
    
    %% Check condition for convergence on the residual vector
    residualNorm = norm(RHS(bc.freeDOF));
    msgNR = sprintf('\t \t ||Fres|| = %d at Newton iteration No. = %d \n',residualNorm,i_newton);
    fprintf([reverseStrNR, msgNR]);
    reverseStrNR = repmat(sprintf('\b'), 1, length(msgNR));
    if residualNorm<=physics.nonLinearScheme.eps
        msgANR = sprintf('\n \t Newton iterations converged! \n \n \n');
        fprintf(msgANR);
        
        % Exit the Newton's loop after convergence
        break;
    end
    
    % If the Newton iterations do not converge after specified limit:
    if i_newton==physics.nonLinearScheme.maxIter
        msgANR = sprintf('\n \t Newton iterations did not converge. \n \n \n');
        fprintf(msgANR);
    end
    
    %% Solve the system for the velocity and the pressure at the new time step
    dup = G(bc.freeDOF, bc.freeDOF) \ RHS(bc.freeDOF);
    up(bc.freeDOF) = up(bc.freeDOF) + dup;
    
    %% Re-assemble to the complete vector of unknowns
    up(bc.drichletDOF) = bc.drichletVector(bc.drichletDOF);
    
% End of Newton Iterations
end

%% Calculating the nodal force vector
forceVector = K*up - F;
fFSI = forceVector(mesh.fsiDOF);


%% Delete and reset the output strings
msgEmpty = sprintf('');
reverseStrANR = repmat(sprintf('\b'), 1, length(msgANR));
reverseStrNR = repmat(sprintf('\b'), 1, length(msgNR));
reverseStrPNR = repmat(sprintf('\b'), 1, length(msgPNR));
fprintf([reverseStrANR, msgEmpty]);
fprintf([reverseStrNR, msgEmpty]);
fprintf([reverseStrPNR, msgEmpty]);


%% Appendix
fprintf('\n');
fprintf('________Nonlinear Transient Analysis Ended________\n');
fprintf('\n');

end