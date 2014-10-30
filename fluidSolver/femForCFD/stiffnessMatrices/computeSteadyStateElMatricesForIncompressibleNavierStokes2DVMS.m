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
function [Ke,KTe,Fe] = computeSteadyStateElMatricesForIncompressibleNavierStokes2DVMS(tau, physics, u, uGp, N, dN, dDotN, laplaceN, Np, dNp)
%% Function documentation
%
% Returns the element stiffness, tangent stiffness and mass matrix as well
% as the element load vector due to body forces using the Variational
% Multiscale stabilization method
%
%% Function main body

%% 3. Compute all necessary matrices and vectors
tauM = tau(1,1); % Stabilization for momentum Equation
tauC = tau(2,2); % Stabilization for Continuity Equation

LopN1 = uGp(1)*dDotN(1)+uGp(2)*dDotN(2);
LopN2 = uGp(1)*dDotN(4)+uGp(2)*dDotN(5);
LopN3 = uGp(1)*dDotN(7)+uGp(2)*dDotN(8);

BoperatorConvection = [ LopN1 0 0 LopN2 0 0 LopN3 0 0;
                        0 LopN1 0 0 LopN2 0 0 LopN3 0;
                        0 0 0 0 0 0 0 0 0 ];
                    
b = physics.bodyForce.magnitude*physics.bodyForce.direction;


%% 4. Compute the diffusion matrix at the Gauss Point
De = physics.nue*(dN'*dN);

% Auxiliary tangent convection B-operator matrix

duN = [ u(1)* dDotN(1)+u(4)*dDotN(4)+u(7)*dDotN(7), u(2)*dDotN(1)+u(5)*dDotN(4)+u(8)*dDotN(7), 0;
        u(1)* dDotN(2)+u(4)*dDotN(5)+u(7)*dDotN(8), u(2)*dDotN(2)+u(5)*dDotN(5)+u(8)*dDotN(8), 0;
        0,                                          0,                                         0 ];
CTangentBOperator = duN*N;

%% 5. Compute the convection matrix at the Gauss Point
Ce = N'*BoperatorConvection;

%% 6. Compute the convection-tangent matrix at the Gauss Point
KTe = (N'*CTangentBOperator);

%% 7. Compute the velocity-pressure coupling matrix at the Gauss Point
Pe = -1*dDotN' * Np;

%% 8. Compute the stabilization matrices on the Gauss Point
%% 8i. Matrices related to diffusion nue*Laplacian(v)

% Diffusion-diffusion stabilization:
Dd = -1*(physics.nue^2)*(laplaceN'*laplaceN);

% Diffusion-convection stabilization:
Dc = physics.nue*(laplaceN'*BoperatorConvection);

% Diffusion-pressure stabilization:
Dp = physics.nue*(laplaceN'*dNp);

%% 8ii. Matrices related to convection v.Gradient(v)

% Convection-diffusion stabilization:
Cd = - Dc';

% Convection-convection stabilization:
Cc = BoperatorConvection'*BoperatorConvection;

% Convection-pressure stabilization:
Cp = BoperatorConvection'*dNp;

%% 8iii. Matrices related to pressure gradient Nabla(p)

% Pressure-diffusion stabilization: (Continuity equation)
Pd = - Dp';

% Pressure-convection stabilization: (Continuity equation)
Pc = Cp';

% Pressure-pressure stabilization: (Continuity equation)
Pp = dNp'*dNp;

%% 8iv. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
Div = dDotN'*dDotN; % (Continuity equation)

%% 9. Compute the mass matrix
Me = N'*N;

%% 10. Compute the stabilization matrices on the mass
%% 10i. Mass-diffusion stabilization
Md = physics.nue*laplaceN'*N;

%% 10ii. Mass-convection stabilization
Mc = BoperatorConvection'*N;

%% 10iii. Mass-pressure stabilization
Mp = dNp'*N;

%% 11. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
    % Stabilization vector related to nue*Laplacian(v) (reaction):
    Fd = physics.nue*laplaceN';
    
    % Stabilization due to pressure gradient Nabla(p): (Continuity equation)
    Fp = dNp';
    
    % 11ii Convection-Body force stabilization
    Fc = -1*(BoperatorConvection');


%% 12. Compute the element stiffness matrix
    % Linear stiffness matrix needed for the residual computation
    Ke = De + Ce + Pe - Pe' + tauM*(Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp) + tauC*Div;
    %Ke = De + Pe - Pe' + tauM*(Dd + Dp + Pd + Pp) + tauC*Div;
%% 14. Compute the stabilized force vector due to source at the Gauss point

% Compute the complete right-hand side vector at the Gauss Point

    % Compute the force vector at the Gauss Point
    Fe = (N' + tauM*(Fd + Fp + Fc))*b;
    
 % End of function
end
