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
function [Ke,Fe] = computeElementMatricesForIncompressibleStokes2DVMS(tau, physics, u, uGp, N, dN, dDotN, laplaceN, Np, dNp)
%% Function documentation
%
% Returns the element stiffness, the element load vector due to body forces using the Variational
% Multiscale stabilization method
%% Function main body

%% 3. Compute all necessary matrices and vectors
tauM = tau(1,1); % Stabilization for momentum Equation
tauC = tau(2,2); % Stabilization for Continuity Equation

%% Obtaining the body force vector
b = physics.bodyForce.magnitude*physics.bodyForce.direction;

%% Computing the Diffusion matrix
D = physics.nue*(dN'*dN);

%% Computing the pressure matrix
P = -dDotN'*Np;

%% Computing the bodyforce Vector
B = N'*b;

%% Computing the coupling matrix between pressure and velocity
Conv = Np'*dDotN;


%% STABILIZATION TERMS

%%  Diffusion-Diffusion Stabilization terms
Dd = -physics.nue*physics.nue*(laplaceN'*laplaceN);

%% Diffusion-Pressure Stabilization
Dp = physics.nue*(laplaceN'*dNp);

%% Pressure-Diffusion Stabilization
Pd = -physics.nue*(dNp'*laplaceN);

%% Pressure-Pressure Stabilization
Pp = dNp'*dNp;

%% Bodyforce-Diffusion Stabilization
Bd = physics.nue*(laplaceN'*b);

%% Bodyforce-Pressure Stabilization
Bp  = dNp'*b;

%% Stabilization of coupling terms
Sc = dDotN'*dDotN;

%% Calculating the Element stiffness matrix
Kactual = D + P + Conv;
Kstabl = tauM*(Dd + Dp + Pd + Pp) + tauC*(Sc);
Ke = Kactual + Kstabl;

%% Calculating the Element force matrix
Factual = B;
Fstabl = tauM*(Bd + Bp);
Fe = Factual + Fstabl;






% End of function
end
