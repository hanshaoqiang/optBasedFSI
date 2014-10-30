%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [fFSI] = formulateFSIforces(pFSI,mesh)
% This function calculates the forces on the FSI interface given as node
% numbers in the mesh
fFSI = zeros(2*length(mesh.fsiNodes),1);
% Loop over all the FSI edges
for i = 1:length(mesh.FSIedge)
    node1 = mesh.FSIedge(i,1);
    node2 = mesh.FSIedge(i,2);
    % find the dx and dy for calculating edge length
    dx = mesh.nodes(node1,1) - mesh.nodes(node2,1);
    dy = mesh.nodes(node1,2) - mesh.nodes(node2,2);
    % find the edge length
    edgeLength = sqrt(dx*dx + dy*dy);
    % Finding the position of node1 and node2 in the fsi nodes list
    fsiNode1 = find(mesh.fsiNodes == node1);
    fsiNode2 = find(mesh.fsiNodes == node2);
    
%     % Force on First Node
%     fFSI(2*fsiNode1-1) = fFSI(2*fsiNode1-1) + edgeLength*pFSI(fsiNode1)*mesh.fsiNormal(fsiNode1,1);     % X Component
%     fFSI(2*fsiNode1) = fFSI(2*fsiNode1) + edgeLength*pFSI(fsiNode1)*mesh.fsiNormal(fsiNode1,2);         % Y Component
%     % Force on Second Node
%     fFSI(2*fsiNode2-1) = fFSI(2*fsiNode2-1) + edgeLength*pFSI(fsiNode2)*mesh.fsiNormal(fsiNode2,1);     % X Component
%     fFSI(2*fsiNode2) = fFSI(2*fsiNode2) + edgeLength*pFSI(fsiNode2)*mesh.fsiNormal(fsiNode2,2);         % Y Component
    
    
    % Force on First Node
    fFSI(2*fsiNode1-1) = fFSI(2*fsiNode1-1) + edgeLength*pFSI(fsiNode1)*-dy;     % X Component
    fFSI(2*fsiNode1) = fFSI(2*fsiNode1) + edgeLength*pFSI(fsiNode1)*dx;         % Y Component
    % Force on Second Node
    fFSI(2*fsiNode2-1) = fFSI(2*fsiNode2-1) + edgeLength*pFSI(fsiNode2)*-dy;     % X Component
    fFSI(2*fsiNode2) = fFSI(2*fsiNode2) + edgeLength*pFSI(fsiNode2)*dx;         % Y Component
    
end
