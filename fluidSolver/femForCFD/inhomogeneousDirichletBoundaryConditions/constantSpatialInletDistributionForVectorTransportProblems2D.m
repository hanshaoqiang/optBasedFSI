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
function scalar = constantSpatialInletDistributionForVectorTransportProblems2D(x,y,z,t)
%% Function documentation
%
% Returns the concetration of the vectorial quantity at the given location 
% y on a vertical (y-directed) line which is defined by a lower and an 
% upper value w.r.t. its y-coordinate. The scalar concetration distribution 
% along this line is assumed to be constant with respect to the space, but 
% varying with respect to the time.
%
%   Input :
%   x,y,z : Cartesian coordinates of the node where the prescribed 
%           (inhomogeneous) Dirichlet boundary conditions are applied
%       t : The time instance where to compute the prescribed Dirichlet 
%           boundary conditions
%
%  Output :
%  scalar : The concetration of the scalar quantity at the given location    
%
%% Function main body

% Maximum concetration of the scalar in the center of the application line
Smax = 40;

% Height of the channel
Height = 2;

% Y-component of the lower edge
ylo = 0;

% Y-component of the upper edge
yup = Height;

% Start time of the simulation
Tstart = 0;

% end time of the simulation
Tend = 1;

% Frequency
omega = 1;

t = 0;

% Compute the scalar value at the given location
if t<=5.13*10^-1
    if y>=ylo && y<=yup
        scalar = Smax*abs(cos(omega*pi*t/(Tend-Tstart)))^2;
    else
        scalar = Smax*abs(cos(omega*pi*t/(Tend-Tstart)))^2;
    end
else
    if y>=ylo && y<=yup
        scalar = Smax*abs(cos(omega*pi*t/(Tend-Tstart)))^2;
    else
        scalar = Smax*abs(cos(omega*pi*t/(Tend-Tstart)))^2;
    end
end

end

