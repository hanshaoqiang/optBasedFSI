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
%   Aditya Ghantasala(M.Sc)            (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [empireMesh] = obtainMeshToSendToEMPIRE(mesh)
%% Function documentation
%
% This function returns a mesh which is suitable to send to EMPIRE.
%INPUT:
%	mesh 		: mesh from which the empire mesh should be formulated.
%
%OUTPUT
%	empireMesh 	: mesh containing entities suitable for EMPIRE.

%% Formulating the necessary variables
empireMesh.numNodes = length(mesh.fsiNodes); % Number of nodes on FSI interface
% Formulating the elements on the FSI interface.
mesh = formulateFSIedges(mesh);
empireMesh.numElems = length(mesh.FSIedge); % Nummber of elements on the FSI interface
% Nodal coordinates of the nodes on the FSI interface
% Format :  
%           |x1 x2 x3 ... | 
%   nodes = |y1 y2 x3 ... |
%           |z1 z2 z3 ... |
empireMesh.nodes            = mesh.nodes(mesh.fsiNodes,:)'; 
empireMesh.nodeIDs          = mesh.fsiNodes;
empireMesh.numNodesPerElem  = 2*ones(length(mesh.FSIedge),1); % Number of elements per elements = 2 for line elements.
empireMesh.elems            = mesh.FSIedge; 

% End of the function
end
