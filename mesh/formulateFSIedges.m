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
function [mesh_s] = formulateFSIedges(mesh)
% This function formulates the edges on the FSI interface and the
% corresponding normals.

mesh.fsiNormal = zeros(length(mesh.fsiNodes),2);
% Looping over all the FSI nodes to find the edges
edgeCount = 1;
for i=1:length(mesh.fsiNodes)
    for j=1:length(mesh.fsiNodes)
        noElemntsPerEdge = 0;
        % Loop over all the elements
        for elem = 1:length(mesh.elements)
            if(~isempty(find(mesh.elements(elem,:) == mesh.fsiNodes(i))) && ~isempty(find(mesh.elements(elem,:) == mesh.fsiNodes(j))) && i ~= j && i>j )
                noElemntsPerEdge = noElemntsPerEdge + 1;
            end
        end
        % if the edge is valid then we  put it in the FSI edges list
        if(noElemntsPerEdge ~= 0 && noElemntsPerEdge <= 2)
            mesh.FSIedge(edgeCount,1) = mesh.fsiNodes(i);
            mesh.FSIedge(edgeCount,2) = mesh.fsiNodes(j);
% %             dx = mesh.nodes(mesh.fsiNodes(i),1) - mesh.nodes(mesh.fsiNodes(j),1);
% %             dy = mesh.nodes(mesh.fsiNodes(i),2) - mesh.nodes(mesh.fsiNodes(j),2);
% %             % Calculating normals
% %             mesh.fsiNormal(i,1) = mesh.fsiNormal(i) + 0.5*dx;
% %             mesh.fsiNormal(i,2) = mesh.fsiNormal(i) + 0.5*dy;
% %             
% %             mesh.fsiNormal((j),1) = mesh.fsiNormal(j) + 0.5*dx;
% %             mesh.fsiNormal((j),2) = mesh.fsiNormal(j) + 0.5*dy;
% %            
            edgeCount =  edgeCount +1;
        end
    end
end

mesh_s = mesh;
% End of the function
end
