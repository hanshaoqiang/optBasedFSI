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
function CPd = computeDisplacedControlPointsIGAVectorTransport2D(CP,scalar,component)
%% Function documentation
%
% Returns the new Control Points with displaced z-component corresponding 
% to the postprocessing of the 2D vector transport problems 
% (convection-diffusion-reaction) 
%
%       Input : 
%          CP : control point coordinates and weights
%      scalar : the scalar concetration field
%   component : 1 : ux (scalar concetration in x-direction)
%               2 : uy (scalar concetration in y-direction)
%               3 : norm(u) (scalar concetration magnitude)
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
            CPd(i,j,3) = CP(i,j,3) + scalar(k);
        elseif component == 2
            CPd(i,j,3) = CP(i,j,3) + scalar(k+1);
        elseif component == 3
            CPd(i,j,3) = sqrt(scalar(k)^2 + scalar(k+1)^2);
        end
        
        % Update counter 
        k = k + 2;
    end
end

end