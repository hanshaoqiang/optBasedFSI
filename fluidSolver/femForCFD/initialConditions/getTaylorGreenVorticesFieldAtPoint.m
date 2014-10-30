%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger   (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                 %
%   Dipl.-Math. Andreas Apostolatos     (andreas.apostolatos@tum.de)       %
%   Aditya Ghantasala                   (aditya.ghantasala@tum.de)        %
%   _______________________________________________________________       %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uvp = getTaylorGreenVorticesFieldAtPoint(physics, geometry, mesh, time, nodeForComparission)
%% Documentation of function
%  This method returns the field for Taylor Green Vortices depending on the
%  time of simulation.
%
%   Input :
%       physics : structure containing the physics of the simulation
%       geometry: structure containing the geometry entities like NURBS
%                   parameters etc
%       mesh    : structure containing the mesh properties like the nodes
%                   and elements.
%       time    : time at which the field is required.
%   Output :
%       uvp     : velocity and pressure field for taylor green vortices.
%% Function main body
uvp = zeros(3,1);
% Obtaining all the x and y coordinates of all the nodes on the mesh. 
x = mesh.nodes(nodeForComparission,1);
y = mesh.nodes(nodeForComparission,2);
f = exp(-2*time*(physics.nue));

%% Calculating the U component
uvp(1:3:end) = (sin(x).*cos(y))*f;

%% Calculating the V component
uvp(2:3:end) = -1*(cos(x).*sin(y))*f;

%% Calculating Pressure
uvp(3:3:end) = 1*0.25*(f^2)*(cos(2*x) + cos(2*y));



% end of the function
end

