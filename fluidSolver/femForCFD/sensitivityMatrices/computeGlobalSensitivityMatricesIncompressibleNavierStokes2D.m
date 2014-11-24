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
%   Reza Najian (M.Sc)                 (reza.najian-asl@tum.de)           %
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de          %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dKdx,dKTangentdx,dMdx,dFdx] = computeGlobalSensitivityMatricesIncompressibleNavierStokes2D(up, mesh, meshVel, physics, transient)
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
dKdx = zeros(ndof,2*length(mesh.nodes));
% Initialize global tangent stiffness matrix
dKTangentdx = zeros(ndof,2*length(mesh.nodes));

% Initialize the global mass matrix
dMdx = zeros(ndof,2*length(mesh.nodes));

% Initializing the body force vector
dFdx = zeros(ndof,2*length(mesh.nodes));

% Obtaining Gauss Intergration points (in parametric space) and weights
[~, gWeight, gCord] = getGaussCordinatesOnTriangle(physics.int.nGauss);


%% 3. loop over all elements
for i = 1:length(mesh.elements)
    %% 3i. Initialize element matrices and force vectors
    
    % Initialize variation for element stiffness matrix
    dKedx = zeros(ndof_e,6);
    
    % Initialize variation for the element tangent stiffness matrix
    dKTedx = zeros(ndof_e,6);
    
    % Initialize variation for the element mass matrix
    dMedx = zeros(ndof_e,6);
   
    % Initialize variation for the element body force vector
    dFedx = zeros(ndof_e,6);
    
    %% 3ii. Initialize the element area size
    
    element_vertices = mesh.elements(i,:);
    element_nodes = [   mesh.nodes(element_vertices(1),:);
                        mesh.nodes(element_vertices(2),:);
                        mesh.nodes(element_vertices(3),:) ];
    
    [h, area] = getTriangularElementSizeAndArea(element_nodes);
    dhdX = getVariationOfSideLengthWithCoord(element_nodes);
    
    %% 3iii. Create an element freedom table
    
    % Element freedom table
    EFT1 = zeros(1,ndof_e); % This is same as original EFT
    EFT2 = zeros(1,6);
    % Assign the entries of the element freedom table
    for k=1:3
        % For Rows
        
        % Velocity DOF
        EFT1(1,3*k-2)    = 3*element_vertices(k)-2;
        EFT1(1,3*k-1)    = 3*element_vertices(k)-1;
        % Pressure DOF
        EFT1(1,3*k)      = 3*element_vertices(k);
        
        % For Columns
        EFT2(1,2*k-1)   = 2*element_vertices(k) -1;
        EFT2(1,2*k)     = 2*element_vertices(k);
    end
    
    %% 3iv. Get the element discrete solution vector of the previous Newton iteration step
    upe = up(EFT1);
    mVele = meshVel(EFT1);
    relVele = upe - mVele;
    vel = [upe(1);upe(2);upe(4);upe(5);upe(7);upe(8)];
    mVel = [mVele(1);mVele(2);mVele(4);mVele(5);mVele(7);mVele(8)];
    %% 3v. Loop over all Quadrature Points for the integration of the element stiffness
    for j = 1:physics.int.nGauss
        
        % Calculate shape functions and their derivatives at the current
        % Gauss point
        [N, dN, dDotN, laplaceN] = getDerivativesOfVelocityShapeFunctions(gCord(:,:,j), element_nodes);
        [Np, dNp] = getDerivativesOfPressureShapeFunctions(gCord(:,:,j), element_nodes);
        
        
        %% 3v.7. Compute the stabilization parameters for the momentum and the continuity equations
        % Calculating velocity at the current Gauss point
        uGp = N*upe;
        mVelGP = N*mVele;
        
        relVelGP = uGp-mVelGP;
        
        % Stabilization for the momentum equation
        if strcmp(physics.stabilization.parametersType,'automatic')
            tauM = ( physics.stabilization.Ct/(transient.dt) + 2*norm(relVelGP)/h + 4*physics.nue/h^2 )^(-1);
%             tauM = ( (physics.stabilization.Ct/transient.dt)^2 + (4*physics.nue/(h^2))^2 + (2*norm(relVelGP)/h)^2 )^(-0.5);  % Based on INRIA paper
            tauM = tauM*physics.stabilization.scalingTauM; % Scalling
            dtauMdX = physics.stabilization.scalingTauM * (-(tauM^2)) * ( 0 + ((-2*norm(relVelGP))/(h^2))*(dhdX) + (-8*physics.nue/(h^3))*dhdX );
            
            % Stabilization for the continuity equation
            tauC = (physics.nue + 0.5*h*norm(relVelGP));
%             tauC = ((physics.nue)^2 + (norm(relVelGP)*h/2)^2)^0.5; % Based on the INRIA Paper
            tauC = tauC*physics.stabilization.scalingTauC; % Scalling
            dtauCdX = physics.stabilization.scalingTauC * (0 + 0.5*norm(relVelGP)*dhdX);
            
            
            tau = [ tauM,0;
                0, tauC];
            
        elseif strcmp(physics.stabilization.parametersType,'manual')
            tauM = physics.stabilization.tauM*physics.stabilization.scalingTauM;
            tauC = physics.stabilization.tauC*physics.stabilization.scalingTauC;
            tau = [ tauM,0;
                0, tauC];
        end
        
        %% 3v.8. Compute the variation of element stiffness, tangent stiffness, mass matrices and body force vector at the quadrature point
        
        % Compute the variation of element stiffness, tangent stiffness,
        % mass matrices and body force vector
        if strcmp(physics.stabilization.type,'variationalMultiscale')
            [dKedxOnGP, dKTedxOnGP, dMedxOnGP, dFedxOnGP] = computeElSensitivities(element_nodes,tau, physics, relVele, relVelGP, N, dN, dDotN, laplaceN, Np, dNp, dtauMdX, dtauCdX, transient);
        end
        
        
        %% Weighing the sensitivities
        dKedx  = dKedx     +   dKedxOnGP*gWeight(j);
        dKTedx = dKTedx    +   dKTedxOnGP*gWeight(j);
        dMedx  = dMedx     +   dMedxOnGP*gWeight(j);
        dFedx  = dFedx     +   dFedxOnGP*gWeight(j);
        
    end
    
    %% Accumulating the sensitivities
    dKdx(EFT1, EFT2)        =           dKdx(EFT1, EFT2)  + area*(dKedx);
    dKTangentdx(EFT1, EFT2) =    dKTangentdx(EFT1, EFT2)  + area*(dKTedx);
    dMdx(EFT1, EFT2)        =           dMdx(EFT1, EFT2)  + area*(dMedx);
    dFdx(EFT1, EFT2)        =           dFdx(EFT1, EFT2)  + area*(dFedx);
end

dKdx        = sparse(dKdx);
dKTangentdx = sparse(dKTangentdx);
dMdx        = sparse(dMdx);
dFdx        = sparse(dFdx);

% End of the function
end

