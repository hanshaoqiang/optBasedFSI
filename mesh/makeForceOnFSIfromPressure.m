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
function [fFSI] = makeForceOnFSIfromPressure(pressureFSI, mesh_s)
% This method calculates forces ont the FSI interface based on the
% integration of pressure on the FSI surface

% This method makes the skewed forces on the FSI interface normal to the
% FSI interface defined as edges in the mesh
fFSI = zeros(2*length(mesh_s.fsiNodes),1);
normals = zeros(length(mesh_s.fsiNodes),2);
unitNormals = zeros(length(mesh_s.fsiNodes),2);
forceMagFSI = zeros(length(mesh_s.fsiNodes),1);

% Loop over all the edges
for i = 1:length(mesh_s.FSIedge)
    
    node1 = mesh_s.FSIedge(i,1);
    node2 = mesh_s.FSIedge(i,2);
    
    fsiNode1 = find(mesh_s.fsiNodes == node1);
    fsiNode2 = find(mesh_s.fsiNodes == node2);
    
    
    dx = mesh_s.nodes(node1,1) - mesh_s.nodes(node2,1);
    dy = mesh_s.nodes(node1,2) - mesh_s.nodes(node2,2);
    lengthOfEdge = sqrt(dx*dx + dy*dy); % Is nothing but distance between the two nodes
    
    forceMagFSI(fsiNode1) = forceMagFSI(fsiNode1) + 0.5*lengthOfEdge*pressureFSI(fsiNode1);
    forceMagFSI(fsiNode2) = forceMagFSI(fsiNode2) + 0.5*lengthOfEdge*pressureFSI(fsiNode2);
    
    normals(fsiNode1,1) = normals(fsiNode1,1) + 0.5*-dy;
    normals(fsiNode2,2) = normals(fsiNode2,2) + 0.5*dx;  
end

% Normalizing the normals calculated above.

for i = 1:length(normals)   
    normal = normals(i,:);
    unitNormal = normal/norm(normal);
    unitNormals(i,:) = unitNormal; 
end


% Making the force normal by doing a dot product and then multiplying with 
for i = 1:length(fFSI)/2
    fFSI(i*2-1) = forceMagFSI(i)*unitNormals(i,1); % X component
    fFSI(i*2)   = forceMagFSI(i)*unitNormals(i,2); % Y component
end



% End of functions
end
