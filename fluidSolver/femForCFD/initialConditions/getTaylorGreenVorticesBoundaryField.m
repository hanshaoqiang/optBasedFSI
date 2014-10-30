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
function uvp = getTaylorGreenVorticesBoundaryField(physics, geometry, mesh, time)
%% Documentation of function
%  This method returns the field for Taylor Green Vortices only on the boundary depending on the
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
uvp = zeros(3*length(mesh.nodesOnBoundary),1);
% Obtaining all the x and y coordinates of all the nodes on the mesh. 
% x = mesh.nodes(mesh.nodesOnBoundary,1);
% y = mesh.nodes(mesh.nodesOnBoundary,2);
f = exp(-2*time*(physics.nue));
% 
% %% Calculating the U component
% uvp(1:3:end) = (sin(x).*cos(y))*f;
% 
% %% Calculating the V component
% uvp(2:3:end) = -1*(cos(x).*sin(y))*f;
% 
% %% Calculating Pressure
% uvp(3:3:end) = -1*0.25*(f^2)*(cos(2.*x) + cos(2.*y));


for n = 1:length(mesh.nodesOnBoundary)
        %uvp_e = getTaylorGreenVorticesFieldAtPoint(physics, geometry, mesh, t, mesh.nodesOnBoundary(n));
        x_n = mesh.nodes(mesh.nodesOnBoundary(n),1);
        y_n = mesh.nodes(mesh.nodesOnBoundary(n),2);
        uvp(3*mesh.nodesOnBoundary(n)-2, 1) = sin(x_n)*cos(y_n)*f;
        uvp(3*mesh.nodesOnBoundary(n)-1, 1) = -1*(cos(x_n)*sin(y_n))*f;
        uvp(3*mesh.nodesOnBoundary(n), 1) = 1*0.25*(f^2)*(cos(2*x_n) + cos(2*y_n));
end


% end of the function
end

