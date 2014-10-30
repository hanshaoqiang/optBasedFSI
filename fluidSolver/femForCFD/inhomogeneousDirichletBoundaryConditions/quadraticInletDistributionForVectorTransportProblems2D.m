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
function scalar = quadraticInletDistributionForVectorTransportProblems2D(x,y,z,t)
%% Function documentation
%
% Returns the concetration of the vectorial quantity at the given location 
% y on a vertical (y-directed) line which is defined by a lower and an 
% upper value w.r.t. its y-coordinate. The scalar concetration distribution 
% along this line is assumed to be quadratic reaching its maximum value at 
% the middle point of this line.
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
Smax = 40.0;

% Height of the channel
Height = 2.0;

% Y-component of the lower edge
% ylo = -Height/2;
ylo = 0;

% Y-component of the upper edge
% yup = Height/2;
yup = Height;

% Start time of the simulation
Tstart = 10;

% end time of the simulation
Tend = 0;

% Fix the time
t = 0;

% Compute the scalar value at the given location
if t<=5.13*10^-1
    if y>=ylo && y<=yup
        scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2*cos(5*pi*t/(Tend-Tstart));
    else
        scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2*cos(5*pi*t/(Tend-Tstart));
    end
else
    if y>=ylo && y<=yup
        scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2*cos(5*pi*t/(Tend-Tstart));
    else
        scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2*cos(5*pi*t/(Tend-Tstart));
    end
end

end

