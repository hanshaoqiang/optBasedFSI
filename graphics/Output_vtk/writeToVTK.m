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
%% Script documentation
%
% Task : Takes the given mesh and produes a VTK output file to be post
% processed with paraview
%   INPUT
%       mesh    : mesh containing then nodal coordinates and etc
%       up      : velocity and pressure (ux, uy, p)T for every node
%   fileName    : the name of the file containing the extension ".vtk"
%

function writeToVTK(mesh,up,fileName)


% Reformulating mesh to be used with the downloaded functions
XY = zeros(length(mesh.nodes),2);
% Extracting x and y coordinates
XY(:,1:2) = mesh.nodes(:,1:2);

% Extracting velocity and pressure values
UVP = zeros(length(mesh.nodes),3);

for i = 1:length(mesh.nodes)
   UVP(i,1) = up(3*i-2);
   UVP(i,2) = up(3*i-1);
   UVP(i,3) = up(3*i);
end 

% Extracting the element connectivity table
eCon = zeros(length(mesh.elements),4);
eCon(:,1:3) = mesh.elements(:,1:3);



%% Writing the file using the downloaded methods

twod_to_vtk ( XY, eCon, UVP, fileName, 'FSI_Simulation');


end