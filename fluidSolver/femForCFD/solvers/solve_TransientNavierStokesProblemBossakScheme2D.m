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
function [up, upRate, fFSI] = solve_TransientNavierStokesProblemBossakScheme2D(mesh, bc, physics, transient, upPrevious, upRatePrevious)
%% Function documentation
%
% Solve the transient Navier-Stokes problem in 2D using the Bossak time
% integration scheme for the temporal discretization
%
%% Function main body

reverseStrNR = '';

nDoFs = 3*length(mesh.nodes);
% Initialize the global flux vector
F = zeros(nDoFs,1);

% Initializing the solution vector
up = upPrevious;

% Initializing the upRate 
upRate = upRatePrevious;

% The time step for the transient analysis
dt = transient.dt;

%% 3vii. Loop over all the Newton-Rapson iterations
for i_newton = 1:physics.nonLinearScheme.maxIter
    %% Compute the tangent stiffness, the stiffness and the mass matrix as well as the residual of the finite element system
    
    % Compute the stiffness, tangent stiffness, mass matrix and body
    % force vector at the current nonlinear iteration step
    [K,KTangent,M,F] = computeMatricesForTransientIncompressibleNavierStokes2D(upPrevious, mesh, mesh.velocity, physics, transient);
    
    % Compute the left-hand side nonlinear stiffness matrix for the
    % correction step
    preFactor = ((1-transient.alphaBeta)/transient.gamma/dt);
    G = preFactor*M + K + KTangent;
    
    % Compute the residual of the transient Navier-Stokes problem at
    % the current nonlinear iteration step
    r = (preFactor*M + K)*up - F - ((1-transient.alphaBeta)/transient.gamma/dt)*M*upPrevious - ...
        ((1-transient.alphaBeta)/transient.gamma - 1)*M*upRate;
    
    %% Updating the right hand side vector to account for drichlet boundary conditons
    
    % Reduce the system matrices and load vectors with respect to all
    % prescribed DoFs. Additionally re-number the inhomogeneous degrees of
    % freedom corresponding to the reduced system
    RHS = - r;
    
    % Compute the generalized right-hand side load vector
    if (i_newton==1)
        % Compute the right-hand side load vector due to the inhomogeneous boundary
        % conditions
        % Update the right-hand side vector using the prescribed values
        RHS = RHS - G*bc.drichletVector;
    end
    
    
    %% Solve the system for the velocity and the pressure at the new time step
    dup = G(bc.freeDOF, bc.freeDOF) \ RHS(bc.freeDOF);
    up(bc.freeDOF) = up(bc.freeDOF) + dup;
    
    % Imposing the boundary conditions on solution vector
    up(bc.drichletDOF) = bc.drichletVector(bc.drichletDOF);
    
    
     %% Check condition for convergence on the residual vector
    residualNorm = norm(RHS(bc.freeDOF));
%     residualNorm = norm(dup);
    if residualNorm<=physics.nonLinearScheme.eps
        % Exit the Newton's loop after convergence
        break;
    end
        
    
end % End of Newton Iterations

  if(residualNorm<=physics.nonLinearScheme.eps)
    msgNR = sprintf('\t \t Newton iterations of Fluid Solver CONVERGED after %d iterations with residual ||Fres|| = %d \n',i_newton,residualNorm);
    fprintf([reverseStrNR, msgNR]);
    reverseStrNR = repmat(sprintf('\b'), 1, length(msgNR));
  elseif(i_newton==physics.nonLinearScheme.maxIter)
    msgNR = sprintf('\t \t Newton iterations of Fluid Solver DID NOT CONVERGE after %d iterations with residual ||Fres|| = %d \n',i_newton,residualNorm);
    fprintf([reverseStrNR, msgNR]);
    reverseStrNR = repmat(sprintf('\b'), 1, length(msgNR));
  
  end


%% Delete and reset the output strings
msgEmpty = sprintf('');
reverseStrNR = repmat(sprintf('\b'), 1, length(msgNR));
fprintf([reverseStrNR, msgEmpty]);


%% Update the acceleration field accordingly
upRate = (1/transient.gamma/dt)*(up - upPrevious) - ((1-transient.gamma)/transient.gamma)*upRatePrevious;


% Because the NS equations are scaled by Density. To find actuall forces
% we scale the pressure back with Density.
forceUn = physics.rho.*(M*upRate + K*up);
fFSI = -1*forceUn(mesh.fsiDOF);


end % End of time loop