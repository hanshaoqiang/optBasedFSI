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
function CPd = computeDisplacedControlPointsIGAIncompressibleFlow2D(CP,nodalUnknowns,component)
%% Function documentation
%
% Returns the new Control Points with displaced z-component corresponding 
% to the postprocessing of a 2D incompressible flow problem
%
%       Input : 
%            CP : control point coordinates and weights
% nodalUnknowns : The vector of the nodal unknowns (velocities + pressures)
%     component : 1 : ux (velocity in x-direction)
%                 2 : uy (velocity in y-direction)
%                 3 : p pressure field
%                 4 : norm(u) (velocity magnitude)
%
%      Output :
%         CPd : Displaced control point coordinates and weights 
%
%% Function main body

% Number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Initialize array of deformed control point location
CPd = CP;

% Initialize counter
k = 1;

% Add the displacements to the array CPd
for j = 1:nv
    for i = 1:nu
        % Displace the z-coordinate of the Control Points accordingly
        if component == 1
            CPd(i,j,3) = CP(i,j,3) + nodalUnknowns(k);
        elseif component == 2
            CPd(i,j,3) = CP(i,j,3) + nodalUnknowns(k+1);
        elseif component == 3
            CPd(i,j,3) = CP(i,j,3) + nodalUnknowns(k+2);
        elseif component == 4
            CPd(i,j,3) = sqrt(nodalUnknowns(k)^2 + nodalUnknowns(k+1)^2);
        end
        
        % Update counter 
        k = k + 3;
    end
end

end