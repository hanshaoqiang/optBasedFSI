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
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function upRate = computeInitialRateConditionForTaylorGreenVorticesFlow2D(physics, geometry, mesh, time)
%% Function documentation
%
% Returns the vector with the initial rate conditions for the case of the 
% Taylor-Green vortices benchmark. Note that for the isogeometric setting 
% only the bilinear elements can be used sue to the non-interpolatory 
% nature of the basis functions for the high-order elements
%
%   Input :
%      CP : The set of control point coordinates and weights.
%
%  Output :
%  upRate : The vector of DoFs containing the initial conditions for the
%           fields
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Degrees of freedom
%   
%% Function main body

%% 0. Read input

%% 1. Loop over all the Degrees of freedom
upRate = zeros(3*length(mesh.nodes),1);
% Obtaining all the x and y coordinates of all the nodes on the mesh. 
x = mesh.nodes(:,1);
y = mesh.nodes(:,2);
f = exp(-2*time*(physics.nue));

%% Calculating the U component
upRate(1:3:end) = -2*physics.nue*(sin(x).*cos(y))*f;

%% Calculating the V component
upRate(2:3:end) = 2*physics.nue*(cos(x).*sin(y))*f;

%% Calculating Pressure
upRate(3:3:end) = -physics.nue*(cos(2*x)+cos(2*y))*(f^2);
    
   
end

