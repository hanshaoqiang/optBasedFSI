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
function p = pressureFieldForTaylorGreenVortices2D(x,y,z,t)
%% Function documentation
%
% Returns the pressure field for the Taylor-Green vortex benchmark
%
%   Input :
%   x,y,z : Cartesian coordinates of the node where the prescribed 
%           (inhomogeneous) Dirichlet boundary conditions are applied
%       t : The time instance where to compute the prescribed Dirichlet 
%           boundary conditions
%
%  Output :
%       u : The pressure field    
%
%% Function main body

% The kinematic viscosity
nue = 1e-3;

% Compute the pressure field
p = -.25*(cos(2*x) + cos(2*y))*exp(-4*t*nue);

end

