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
function up = computeInitialConditionForTaylorGreenVorticesFlow2D(CP)
%% Function documentation
%
% Returns the vector with the initial conditions for the case of the Taylor
% -Green vortices benchmark. Note that for the isogeometric setting only
% the bilinear elements can be used sue to the non-interpolatory nature of
% the basis functions for the high-order elements
%
%   Input :
%      CP : The set of control point coordinates and weights.
%
%  Output :
%      up : The vector of DoFs containing the initial conditions for the
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

% Number of Control Points
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
nCPs = nu*nv;

% Number of degees of freedom
nDoFs = 3*nCPs;

% Initialize output array
up = zeros(nDoFs,1);

%% 1. Loop over all the Degrees of freedom
for l=1:nDoFs
    % Get the numbering of the corresponding Control Point
    p = ceil(l/3);
    j = ceil(p/nu);
    i = p - (j-1)*nu;
    dir = l - ((j-1)*nu + i-1)*3;
    
    % Get the Cartesian coordinates of the Control Point
    x = CP(i,j,1);
    y = CP(i,j,2);
    
    % Assign the initial condition according to the nature of the DoF (i.e. 
    % if it velocity in x or y-direction or pressure)
    if dir==1 % velocity in x-direction
        up(l) = -cos(x)*sin(y);
    elseif dir==2 % velocity in y-direction
        up(l) = sin(x)*cos(y);
    elseif dir==3 % pressure
        up(l) = -.25*(cos(2*x) + cos(2*y));
    end
end



end

