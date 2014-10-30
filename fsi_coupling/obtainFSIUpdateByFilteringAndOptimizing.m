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
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method computes the updates for fsi velocity and fsi
% positon/coordinates based on the optimization and filtering approach
%   Input:
%       dI_dX: 				: The variation of objective function with respect
%							  to the design variables.
%
%
%   Output:
%       fsi_cord_update     : The update to be added/applied to the
%                             coordinates of the the FSI interface. 
%       fsi_vel_update      : The update to be added/applied to the
%                             velocity of the FSI interface. This done to
%                             be consistant with all the terms that are
%                             applied on the fluid solver.
%
%
%       The displacement is useful for the calculation of the surrounding
%       velocity distribution, and the velocity is used as a Drichlet
%       boundary condition on the FSI interface. So they both should be
%       consistant in formulation.
%
function [fsi_cord_update, fsi_vel_update, df_dz] = obtainFSIUpdateByFilteringAndOptimizing(mesh_fsi, vel_fsi,...
                                                                                     objectiveFunction_disp, controlField_disp, ...
                                                                                     objectiveFunction_vel, controlField_vel)

% Degree of freedoms on the interface
nDOF = 2*length(mesh_fsi.interfaceNodeNumbers);
% Empty global matrices declaration
Bcord = zeros(nDOF, nDOF);
Bvel = zeros(nDOF, nDOF);
df_dz = zeros(nDOF, 1);
df_dvel = zeros(nDOF, 1);

% Misc properties
nGauss = 2; %%% Number of Gauss points used for integration 2 or 3

% Obtaining Gauss Intergration points (in parametric space) and weights
[gaussWeights, gaussCord] = getGaussCordinatesOn1DElement(nGauss);


%% Computing the element wise Bcord, Bvel, df_dz, df_dvel and forming the 
%%% Global matrix.
%%% like in the classical finite elements we first loop over all the
%%% elements on the fsi interface, and then add the contributions, on the
%%% gauss points.

% Loop over all the elements on the FSI interface
for elem = 1:length(mesh_fsi.elements)
    
    % Element matrices intialization
    Bcord_e = zeros(4,4);
    Bvel_e = zeros(4,4);
    
    % Node numbers making up the current element
    element_vertices = mesh_fsi.elements(elem,:);
    % Obtaining the nodes corresponding to the current element
    element_nodes = [   mesh_fsi.nodes_base(element_vertices(1),:);
                        mesh_fsi.nodes_base(element_vertices(2),:)];
    node1_position = find(mesh_fsi.interfaceNodeNumbers == element_vertices(1));
    node2_position = find(mesh_fsi.interfaceNodeNumbers == element_vertices(2));
    %% CHECK 
    if(node1_position == 4 || node1_position == 3 || node1_position == 2 || node1_position == 1)
    end
    if(node2_position == 4 || node2_position == 3 || node2_position == 2 || node2_position == 1)
    end
                    
    % Objective function for the current element
    objFun_disp_e = [sqrt(abs(objectiveFunction_disp(2*node1_position -1))^2 + abs(objectiveFunction_disp(2*node1_position))^2);
                     sqrt(abs(objectiveFunction_disp(2*node2_position -1))^2 + abs(objectiveFunction_disp(2*node2_position))^2) ];
    
    % Objective function for the current element
    objFun_vel_e = [sqrt(abs(objectiveFunction_vel(2*node1_position -1))^2 + abs(objectiveFunction_vel(2*node1_position))^2);
                     sqrt(abs(objectiveFunction_vel(2*node2_position -1))^2 + abs(objectiveFunction_vel(2*node2_position))^2) ];

    % Control field for the current element
    contField_disp_e = [controlField_disp(2*node1_position-1); controlField_disp(2*node1_position);  controlField_disp(2*node2_position-1); controlField_disp(2*node2_position)];
    contField_vel_e = [controlField_vel(2*node1_position-1); controlField_vel(2*node1_position);  controlField_vel(2*node2_position-1); controlField_vel(2*node2_position)];
    
    % Velocity on the nodes of the current elements.
    vel_e(1,:) = [vel_fsi(2*node1_position-1) vel_fsi(2*node1_position)];
    vel_e(2,:) = [vel_fsi(2*node2_position-1) vel_fsi(2*node2_position)];
    
    %%% Forming the Element freedom table
    EFT = zeros(2*2,1);
    % First node
    EFT(1) = 2*node1_position-1;
    EFT(2) = 2*node1_position;
    % Second node
    EFT(3) = 2*node2_position-1;
    EFT(4) = 2*node2_position;
        
    % Lenghtn ofthe element to be used for converting the integeral from
    % the parametric coordinates to the real coordinates
    x1 = element_nodes(1,1); y1 = element_nodes(1,2);
    x2 = element_nodes(2,1); y2 = element_nodes(2,2);
    
    length_x = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    
    %%% Loop over all the Gaussian points
    for gaussPt = 1:nGauss
        
        % Obtaining the Filter function values at the gaussian point.
        %% It is important that how we formulate the N matrix in the classical FEM.
        %%% It should be nDOF x nDOF for each element.
        %%% Here we also treat for the special Filter functions which are
        %%% required for the nodes on the domain edges. The filter
        %%% functions should have radius that the filter function should
        %%% not go beyond the domain. So we decrease the radius of the
        %%% filter function. The reference is in Bletzinger's paper.
        
        %%% TODO :: Modify the filter function using the radius of filter. 
        %%%         How do we know which nodal filter is spanning the
        %%%         current element. Because if the radius is big the
        %%%         filter function can span multiple elements before it
        %%%         becomes zero. The same case as nurbs basis functions.
        %%%         Observe from Andrea's solver.
        
        Fcord   = obtainGeometryFilterFunctionValueAtGaussPoint(gaussCord(:,:,gaussPt), element_nodes);
        Fvel    = obtainVelocityFilterFunctionValueAtGaussPoint(gaussCord(:,:,gaussPt), element_nodes);
                
        % Obtaining the Control field (here pressure) values at the
        % gaussian point.
        P_disp = obtainShapeFunctionValueAtGaussPoint(gaussCord(:,:,gaussPt),element_nodes, contField_disp_e);
        P_vel = obtainShapeFunctionValueAtGaussPoint(gaussCord(:,:,gaussPt),element_nodes, contField_vel_e);
       
        Bcord_e_gp = Fcord'*P_disp;       %%% 
        Bvel_e_gp  = Fvel'*P_vel;        %%% If we take same filter functions for coordinates 
                                     %%% and velocity, this will be equal to Bcord_e_gp
            
        % Adding the gaussian contribution to the element matrix
        Bcord_e = Bcord_e + Bcord_e_gp * gaussWeights(gaussPt); 
        Bvel_e  = Bvel_e  + Bvel_e_gp  * gaussWeights(gaussPt);
 
    end   
    
    % Obtaining the variation of the objective function with respect
    % to the design variables coordinates and velocities.                
    df_dz_e     = obtainVariationOfObjWithGeometry( element_nodes, objFun_disp_e );  %%%
    df_dvel_e   = obtainVariationOfObjWithVelocity( vel_e, element_nodes,  objFun_vel_e );  %%% 
    
    %%% Obtaining the damping factor for the nodes of current element. 
    %%% The nodes which are away from any of the nodes on the edge will not
    %%% have any effect. But the nodes which are near to the edges will
    %%% have to be damped so they will be effected. We multiply the
    %%% variation of the objective function w.r.t z and velocity with this
    %%% damping factor. The formula for the damping factor is from the
    %%% implementation in CARAT File "ShapeBasis.cpp" line 341.
    
    %                                   pi
    % dampingFactor = 1 - cos( --------------------- * (nodeDistance - fixedDistance) )
    %                           2 * dampingDistance
    
    % fixedDistance     : is the distance from the edge till where the
    %                     edge effect is considered. 
    % dampingDistance   : Is the distance from the edge till where the
    %                     damping is applied.
    % nodeDistance      : is the distance of the current from the nearest
    %                     node on the edge. 
    fixedDistance = 0.05;
    dampingDistance = 0.05;
    dampingFactor = obtainDampingFactorForElementNodes(mesh_fsi, element_nodes, [node1_position, node2_position], fixedDistance, dampingDistance);
    
    % Applying the damping factor to the variations
    df_dz_e     = dampingFactor .* df_dz_e;
    df_dvel_e   = dampingFactor .* df_dvel_e;
       
    
    %%% Assembling the element contribution to the global matrix using EFT
    Bcord(EFT,EFT)  =  Bcord(EFT,EFT)   + Bcord_e   * length_x;
    Bvel(EFT,EFT)   =  Bcord(EFT,EFT)   + Bvel_e    * length_x; 
    df_dz(EFT)      =  df_dz(EFT)       + df_dz_e;   
    df_dvel(EFT)    =  df_dvel(EFT)     + df_dvel_e;
    
    
% End of loop over elements
end

%% Computing the geometry update
%%% Here 
%%% Bcord   : the global (on the fsi interface) B matrix 
%%% df_dz   : the global (on the fsi interface) df_dz
fsi_cord_update = Bcord*eye(length(Bcord))*Bcord' * df_dz;



%% Computing the velocity update
%%% Here 
%%% Bcord   : the global (on the fsi interface) B matrix 
%%% df_dvel : the global (on the fsi interface) df_dvel
%%%           That is the variation of the objective function with respect
%%%           to the design parameter, here since we want to optimize the
%%%           velocity, the design parameter will be velocity.
fsi_vel_update = Bvel*eye(length(Bvel))*Bvel' * df_dvel;



end
