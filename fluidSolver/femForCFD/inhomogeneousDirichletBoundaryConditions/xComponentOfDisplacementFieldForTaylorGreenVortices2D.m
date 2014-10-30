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
function u = xComponentOfDisplacementFieldForTaylorGreenVortices2D(x,y,z,t)
%% Function documentation
%
% Returns the x-component of the velocity field for the Taylor-Green vortex 
% benchmark
%
%   Input :
%   x,y,z : Cartesian coordinates of the node where the prescribed 
%           (inhomogeneous) Dirichlet boundary conditions are applied
%       t : The time instance where to compute the prescribed Dirichlet 
%           boundary conditions
%
%  Output :
%       u : The x-component of the velocity field    
%
%% Function main body

% The kinematic viscosity
nue = 1e-3;

% Compute the x-component of the velocity field
u = -cos(x)*sin(y)*exp(-2*t*nue);

end

