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
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,KTangent,F] = computeSteadyStateMatricesForIncompressibleNavierStokes2D(up, mesh, physics)
%% Function documentation
%
% Returns the stiffness, the tangent stiffness and the mass matrices for
% the incompressible transient Navier-Stokes problem in 2D.
%
%% Function main body

%% 0. Read input

% Total number of degrees of freedom for the 2D vector transport equation
ndof = 3*length(mesh.nodes);

% Number of DoFs affected the element under study
ndof_e = 3*3;

% Initialize global stiffness matrix
K = zeros(ndof,ndof);

% Initialize global tangent stiffness matrix
KTangent = zeros(ndof,ndof);

% Initializing the body force vector
F = zeros(ndof,1);

% Obtaining Gauss Intergration points (in parametric space) and weights
[~, gWeight, gCord] = getGaussCordinatesOnTriangle(physics.int.nGauss);


%% 3. loop over all elements
for i = 1:length(mesh.elements)
    %% 3i. Initialize element matrices and force vectors
    
    % Initialize element stiffness matrix
    Ke = zeros(ndof_e,ndof_e);
    
    % Initialize the element tangent stiffness matrix
    KTe = zeros(ndof_e,ndof_e);
    
    % Initialize the element body force vector
    Fe = zeros(ndof_e,1);
    
    %% 3ii. Initialize the element area size
    
    element_vertices = mesh.elements(i,:);
    element_nodes = [   mesh.nodes(element_vertices(1),:);
        mesh.nodes(element_vertices(2),:);
        mesh.nodes(element_vertices(3),:) ];
    
    [h, area] = getTriangularElementSizeAndArea(element_nodes);
    h = 100*h;
    
    %% 3iii. Create an element freedom table
    
    % Element freedom table
    EFT = zeros(1,ndof_e);
    % Assign the entried of the element freedom table recursively
    for k=1:3
        % Velocity DOF
        EFT(1,3*k-2)    = 3*element_vertices(k)-2;
        EFT(1,3*k-1)    = 3*element_vertices(k)-1;
        EFT(1,3*k)      = 3*element_vertices(k);
    end
    
    %% 3iv. Get the element discrete solution vector of the previous Newton iteration step
    upe = up(EFT);
    vel = [upe(1);upe(2);upe(4);upe(5);upe(7);upe(8)];
    
    %% 3v. Loop over all Quadrature Points for the integration of the element stiffness
    for j = 1:physics.int.nGauss
        
        % Calculate shape functions and their derivatives at the current
        % Gauss point
        [N, dN, dDotN, laplaceN] = getDerivativesOfVelocityShapeFunctions(gCord(:,:,j), element_nodes);
        [Np, dNp] = getDerivativesOfPressureShapeFunctions(gCord(:,:,j), element_nodes);
        
        
        %% 3v.7. Compute the stabilization parameters for the momentum and the continuity equations
        
        % Calculating velocity at the current Gauss point
        uGp = N*upe;
        
        % Stabilization for the momentum equation
        if strcmp(physics.stabilization.parametersType,'automatic')
            tauM = ( (norm(uGp)/h) + (4*physics.nue/(h^2)) )^(-1);
%             tauM = ( (4*physics.nue/(h^2))^2 + (4*norm(uGp)/h)^2 )^(-0.5);
            tauM = tauM*physics.stabilization.scalingTauM; % Scalling
            
            % Stabilization for the continuity equation
            tauC = (physics.nue + 0.5*h*norm(uGp));
%             tauC = ((physics.nue)^2 + (norm(uGp)/h)^2)^0.5;
            tauC = tauC*physics.stabilization.scalingTauC; % Scalling
            
            tau = [ tauM,0;
                0, tauC];
            
        elseif strcmp(stabilization.parametersType,'manual')
            tauM = stabilization.tauM*stabilization.scalingTauM;
            tauC = stabilization.tauC*stabilization.scalingTauC;
            tau = [ tauM,0;
                    0, tauC];
        end
        
        %% 3v.8. Compute the element stiffness, tangent stiffness, mass matrices and body force vector at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
        
        % Compute the element stiffness, tangent stiffness,
        % mass matrices and body force vector
        if strcmp(physics.stabilization.type,'variationalMultiscale')
            [KeOnGP,KTeOnGP,FeOnGP] = computeSteadyStateElMatricesForIncompressibleNavierStokes2DVMS(tau, physics, upe, uGp, N, dN, dDotN, laplaceN, Np, dNp);
        end
        
        %% 3v.9. Add the contribution from the Gauss Point
        Ke = Ke + KeOnGP*gWeight(j);
        KTe = KTe + KTeOnGP*gWeight(j);
        Fe = Fe + FeOnGP*gWeight(j);
    end
    
    
    %% 3vii. insert Ke into K via the element freedom table
    K(EFT, EFT) =           K(EFT, EFT)         + area*(Ke);
    KTangent(EFT, EFT) =    KTangent(EFT, EFT)  + area*(KTe);
    F(EFT) =                F(EFT)              + area*(Fe);
end

K = sparse(K);
KTangent = sparse(KTangent);
F = sparse(F);

% End of the function
end

