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
function flux = sinusoidalOutletDistributionForScalarTransportProblems2D(x,y,z,t)
%% Function documentation
%
% Returns the flux of the scalar quantity at the given location y
% on a vertical (y-directed) line which is defined by a lower and an upper
% value w.r.t. its y-coordinate and at a time instance t \in [Tstart, Tend] 
% The flux of the scalar distribution along this line is assumed to be 
% sinusoidal with respect to both time and space.
%
%  Input :
%  x,y,z : The Cartesian coordinates of the node where to evaluate the flux
%      t : The time instance where the flux to be computed
%
% Output :
%   flux : The flux of the scalar quantity at the given location    
%
%% Function main body

% Maximum concetration of the scalar in the center of the application line
Fmax = 4;

% Height of the channel
Height = 10;

% Y-component of the lower edge
ylo = -Height/2;

% Y-component of the upper edge
yup = Height/2;

% Start time of the simulation
Tstart = 10;

% end time of the simulation
Tend = 0;

% Compute the applied flux value at the given location
flux = -10*Fmax*yup*ylo*sin(pi*y/((yup-ylo))/2)*sin(pi*t/(Tend-Tstart)/4)/(yup-ylo)^2;

end

