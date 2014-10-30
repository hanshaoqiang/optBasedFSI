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
% Obtaining the damping factor for the nodes of current element.
% The nodes which are away from any of the nodes on the edge will not
% have any effect. But the nodes which are near to the edges will
% have to be damped so they will be effected. We multiply the
% variation of the objective function w.r.t z and velocity with this
% damping factor. The formula for the damping factor is from the
% implementation in CARAT File "ShapeBasis.cpp" line 341.
%
%                                   pi
% dampingFactor = 1 - cos( --------------------- * (nodeDistance - fixedDistance) )
%                           2 * dampingDistance
%   Input:
%       mesh_fsi          : The mesh of the FSI interface.
%       element_nodes     : The nodal coordinates of the nodes making up the element
%       node_numbers      : The numbers of the nodes that make up the
%                           element.
%       fixedDistance     : is the distance from the edge till where the
%                           edge effect is considered.
%       dampingDistance   : Is the distance from the edge till where the
%                           damping is applied.
%
%   Output: 
%       dampingFactor     : The vector of the damping factors for all DOFs
%                           of the given element.
%
function [dampingFactor] = obtainDampingFactorForElementNodes(mesh_fsi, element_nodes, node_numbers, fixedDistance, dampingDistance)

% Calculating the distances of the given nodes from the nearest edge node
% (provided in the mesh_fsi)
% This will be the minimum of all the distances from the nodes on the edges
dampingFactor = zeros(4,1);
x1 = element_nodes(1,1); y1 = element_nodes(1,2);
x2 = element_nodes(2,1); y2 = element_nodes(2,2);
for i = 1:length(mesh_fsi.edgeNodes)
    x = mesh_fsi.nodes(mesh_fsi.edgeNodes(i),1);
    y = mesh_fsi.nodes(mesh_fsi.edgeNodes(i),2);
    
    lengthx(i) = sqrt( (x1-x)*(x1-x) + (y1-y)*(y1-y) );
    lengthy(i) = sqrt( (x2-x)*(x2-x) + (y2-y)*(y2-y) );
end

nodeDistance  = [min(lengthx); min(lengthy); min(lengthx); min(lengthy)];
% Calculating the damping factor
for i = 1:length(nodeDistance)
    if(nodeDistance(i) <= fixedDistance)
        dampingFactor(i) = 0;
    elseif(nodeDistance(i) > fixedDistance + dampingDistance)
         dampingFactor(i) = 1;
    else 
        dampingFactor(i) = 1 - cos((pi./(2 * dampingDistance)) .* (nodeDistance(i) - fixedDistance) ) ;
    end
end

% End of the function
end

