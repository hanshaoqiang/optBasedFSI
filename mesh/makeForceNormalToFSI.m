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
function [fFSI] = makeForceNormalToFSI(fFSIin, mesh)
% This method makes the skewed forces on the FSI interface normal to the
% FSI interface defined as edges in the mesh
fFSI = zeros(length(fFSIin),1);
normals = zeros(length(mesh.fsiNodes),2);
unitNormals = zeros(length(mesh.fsiNodes),2);

% Loop over all the edges
for i = 1:length(mesh.FSIedge)
    
    node1 = mesh.FSIedge(i,1);
    node2 = mesh.FSIedge(i,2);
    
    fsiNode1 = find(mesh.fsiNodes == node1);
    fsiNode2 = find(mesh.fsiNodes == node2);
    
    
    dx = mesh.nodes(node1,1) - mesh.nodes(node2,1);
    dy = mesh.nodes(node1,2) - mesh.nodes(node2,2);
    
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
for i = 1:length(fFSIin)/2

    normal = unitNormals(i,:);
    forceVector = [fFSIin(2*i-1),fFSIin(2*i)];
    
    projection = normal(1)*forceVector(1) + normal(2)*forceVector(2);
 
    fFSI(i*2-1) = projection*unitNormals(i,1);
    fFSI(i*2)   = projection*unitNormals(i,2);
end


% End of method
end