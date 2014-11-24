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
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de)         %
%   Reza Najian (M.Sc)                 (reza.najian-asl@tum.de)           %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dKedX, dKTedX, dMedX, dFedX] = computeElSensitivities(nodes, tau, physics, u, uGp, N, dN, dDotN, laplaceN, Np, dNp, dtauMdX, dtauCdX, transient)
%% Function documentation
%
% Returns the element stiffness, tangent stiffness and mass matrix as well
% as the element load vector due to body forces using the Variational
% Multiscale stabilization method
%
%% Function main body

x1 = nodes(1,1);
y1 = nodes(1,2);
% NODE 2
x2 = nodes(2,1);
y2 = nodes(2,2);
% NODE 3
x3 = nodes(3,1);
y3 = nodes(3,2);

% Computing area
area = 0.5*abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );

%% Calculating all the derivatives of the Shape functions
% NOTE : Basis Function used are of Order 1
% #############

% Differential of Basis Function 1
dN1_dx = (y2-y3)/(2*area);
dN1_dy = (x3-x2)/(2*area);

% Differential of Basis Function 2
dN2_dx = (y3-y1)/(2*area);
dN2_dy = (x1-x3)/(2*area);

% Differential of Basis Function 3
dN3_dx = (y1-y2)/(2*area);
dN3_dy = (x2-x1)/(2*area);
dDotN = [dN1_dx dN1_dy 0 dN2_dx dN2_dy 0 dN3_dx dN3_dy 0];
%% Calculatingt the variation of shape funtions with coordinates
%%% This will be a 6x6 matrix
dNdX = getVariationOfVelocityShapeFunctionsWithCord(uGp, nodes);
%%% QUESTION ... There should be variation of pressure shape function
%%% also.Correct ... ??


%% Predefining the variables necessary
dN1dx = dDotN(1);
dN1dy = dDotN(2);
dN2dx = dDotN(4);
dN2dy = dDotN(5);
dN3dx = dDotN(7);
dN3dy = dDotN(8);



%% Fourmulating the matrix dN
% Its a 9x9 matrix
dN = [dN1dx 0 0 dN2dx 0 0 dN3dx 0 0;
      0 dN1dx 0 0 dN2dx 0 0 dN3dx 0;
      0 0 0 0 0 0 0 0 0;
      dN1dy 0 0 dN2dy 0 0 dN3dy 0 0;
      0 dN1dy 0 0 dN2dy 0 0 dN3dy 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 ];
  
dNp = [0 0 dN1dx 0 0 dN2dx 0 0 dN3dx;
       0 0 dN1dy 0 0 dN2dy 0 0 dN3dy;
       0 0 0      0 0 0      0 0 0];     

%
u1x = u(1);
u1y = u(2);
p1  = u(3);
u2x = u(4);
u2y = u(5);
p2  = u(6);
u3x = u(7);
u3y = u(8);
p3  = u(9);

%
Gp1 = Np(3);
Gp2 = Np(6);
Gp3 = Np(9);

% 
ux = uGp(1);
uy = uGp(2);


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


duN = [ u(1)* dDotN(1)+u(4)*dDotN(4)+u(7)*dDotN(7), u(2)*dDotN(1)+u(5)*dDotN(4)+u(8)*dDotN(7), 0;
    u(1)* dDotN(2)+u(4)*dDotN(5)+u(7)*dDotN(8), u(2)*dDotN(2)+u(5)*dDotN(5)+u(8)*dDotN(8), 0;
    0,                                          0,                                         0 ];
CTangentBOperator = duN*N;
% %% FD first 4 matrics
De_b = (physics.nue*(dN'*dN)) * u;
Ce_b = (N'*BoperatorConvection) * u;
Pe_b = (-1*dDotN' * Np) * u;
% val_b = De_b + Ce_b + Pe_b;
Dd_b = (-1*(physics.nue^2)*(laplaceN'*laplaceN))*u;
Dc_b = (physics.nue*(laplaceN'*BoperatorConvection))*u;
Dp_b = (physics.nue*(laplaceN'*dNp))*u;
Cd_b = (- (physics.nue*(laplaceN'*BoperatorConvection))')*u;
Cc_b = (BoperatorConvection'*BoperatorConvection)*u;
Cp_b = (BoperatorConvection'*dNp)*u;
Pd_b = (- (physics.nue*(laplaceN'*dNp))')*u;
Pp_b = (dNp'*dNp)*u;
Pc_b = ((BoperatorConvection'*dNp)')*u;
Div_b = (dDotN'*dDotN)*u;
tauM_b = tauM; tauC_b = tauC;
Mp_b = (dNp'*N)*u;
Mc_b = (BoperatorConvection'*N) * u;

val_b =  De_b + Ce_b + Pe_b + tauM_b*(Dd_b+ Dc_b + Dp_b + Cd_b + Cc_b + Cp_b + Pd_b + Pc_b + Pp_b) + tauC_b*Div_b;
mval_b = tauM_b*(Mp_b + Mc_b);
delta = 1E-8;
% dDotN(4) = dDotN(4) + delta;
% dN2dx = dN2dx + delta;
y1 = y1 + delta;


%% purt recalculation
% Computing area
area = 0.5*abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );


% Differential of Basis Function 1
dN1_dx = (y2-y3)/(2*area);
dN1_dy = (x3-x2)/(2*area);

% Differential of Basis Function 2
dN2_dx = (y3-y1)/(2*area);
dN2_dy = (x1-x3)/(2*area);

% Differential of Basis Function 3
dN3_dx = (y1-y2)/(2*area);
dN3_dy = (x2-x1)/(2*area);
dDotN = [dN1_dx dN1_dy 0 dN2_dx dN2_dy 0 dN3_dx dN3_dy 0];

dN1dx = dDotN(1);
dN1dy = dDotN(2);
dN2dx = dDotN(4);
dN2dy = dDotN(5);
dN3dx = dDotN(7);
dN3dy = dDotN(8);

dN = [dN1dx 0 0 dN2dx 0 0 dN3dx 0 0;
      0 dN1dx 0 0 dN2dx 0 0 dN3dx 0;
      0 0 0 0 0 0 0 0 0;
      dN1dy 0 0 dN2dy 0 0 dN3dy 0 0;
      0 dN1dy 0 0 dN2dy 0 0 dN3dy 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 ];
  
dNp = [0 0 dN1dx 0 0 dN2dx 0 0 dN3dx;
       0 0 dN1dy 0 0 dN2dy 0 0 dN3dy;
       0 0 0      0 0 0      0 0 0]; 
LopN1 = uGp(1)*dDotN(1)+uGp(2)*dDotN(2);
LopN2 = uGp(1)*dDotN(4)+uGp(2)*dDotN(5);
LopN3 = uGp(1)*dDotN(7)+uGp(2)*dDotN(8);

BoperatorConvection = [ LopN1 0 0 LopN2 0 0 LopN3 0 0;
    0 LopN1 0 0 LopN2 0 0 LopN3 0;
    0 0 0 0 0 0 0 0 0 ];
b = physics.bodyForce.magnitude*physics.bodyForce.direction;


duN = [ u(1)* dDotN(1)+u(4)*dDotN(4)+u(7)*dDotN(7), u(2)*dDotN(1)+u(5)*dDotN(4)+u(8)*dDotN(7), 0;
    u(1)* dDotN(2)+u(4)*dDotN(5)+u(7)*dDotN(8), u(2)*dDotN(2)+u(5)*dDotN(5)+u(8)*dDotN(8), 0;
    0,                                          0,                                         0 ];
CTangentBOperator = duN*N;

% EOF purt recalculation


De_a = (physics.nue*(dN'*dN) ) * u;
Ce_a = (N'*BoperatorConvection) * u;
Pe_a = (-1*dDotN' * Np) * u;
% val_a = De_a + Ce_a + Pe_a;
% 
Dd_a = (-1*(physics.nue^2)*(laplaceN'*laplaceN))*u;
Dc_a = (physics.nue*(laplaceN'*BoperatorConvection))*u;
Dp_a = (physics.nue*(laplaceN'*dNp))*u;
Cd_a = (- (physics.nue*(laplaceN'*BoperatorConvection))')*u;
Cc_a = (BoperatorConvection'*BoperatorConvection)*u;
Cp_a = (BoperatorConvection'*dNp)*u;
Pd_a = (- (physics.nue*(laplaceN'*dNp))')*u;
Pp_a = (dNp'*dNp)*u;
Pc_a = ((BoperatorConvection'*dNp)')*u;
Div_a = (dDotN'*dDotN)*u;
Mp_a = (dNp'*N)*u;
Mc_a = (BoperatorConvection'*N) * u;

nodes(1,2) = y1;
[h, area] = getTriangularElementSizeAndArea(nodes);
tauM_a = ( physics.stabilization.Ct/(transient.dt) + 2*norm(uGp)/h + 4*physics.nue/h^2 )^(-1);
tauC_a = (physics.nue + 0.5*h*norm(uGp));
val_a = De_a + Ce_a + Pe_a + tauM_a*(Dd_a+ Dc_a + Dp_a + Cd_a + Cc_a + Cp_a + Pd_a + Pc_a + Pp_a) + tauC_a*Div_a;
mval_a = tauM_a*(Mp_a + Mc_a);

% % val_a = De_a + Ce_a + Pe_a + tauM*(Dd_a+ Dc_a + Dp_a + Cd_a + Cc_a + Cp_a + Pd_a + Pc_a + Pp_a) + tauC*Div_a;
% val_a = tauM_a * (Dd_b+ Dc_b + Dp_b + Cd_b + Cc_b + Cp_b + Pd_b + Pc_b + Pp_b);

FD_taU = (tauM_a-tauM)/delta;

FD = (val_a - val_b)/delta;
FDM = (mval_a - mval_b)/delta;
% FD = (Pp_a - Pp_b)/delta;


%% taking back
% dDotN(4) = dDotN(4) - delta;
% dN2dx = dN2dx - delta;
y1 = y1 - delta;
nodes(1,2) = y1;
area = 0.5*abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );


% Differential of Basis Function 1
dN1_dx = (y2-y3)/(2*area);
dN1_dy = (x3-x2)/(2*area);

% Differential of Basis Function 2
dN2_dx = (y3-y1)/(2*area);
dN2_dy = (x1-x3)/(2*area);

% Differential of Basis Function 3
dN3_dx = (y1-y2)/(2*area);
dN3_dy = (x2-x1)/(2*area);
dDotN = [dN1_dx dN1_dy 0 dN2_dx dN2_dy 0 dN3_dx dN3_dy 0];

dN1dx = dDotN(1);
dN1dy = dDotN(2);
dN2dx = dDotN(4);
dN2dy = dDotN(5);
dN3dx = dDotN(7);
dN3dy = dDotN(8);

dN = [dN1dx 0 0 dN2dx 0 0 dN3dx 0 0;
      0 dN1dx 0 0 dN2dx 0 0 dN3dx 0;
      0 0 0 0 0 0 0 0 0;
      dN1dy 0 0 dN2dy 0 0 dN3dy 0 0;
      0 dN1dy 0 0 dN2dy 0 0 dN3dy 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 ];
  
dNp = [0 0 dN1dx 0 0 dN2dx 0 0 dN3dx;
       0 0 dN1dy 0 0 dN2dy 0 0 dN3dy;
       0 0 0      0 0 0      0 0 0]; 
LopN1 = uGp(1)*dDotN(1)+uGp(2)*dDotN(2);
LopN2 = uGp(1)*dDotN(4)+uGp(2)*dDotN(5);
LopN3 = uGp(1)*dDotN(7)+uGp(2)*dDotN(8);

BoperatorConvection = [ LopN1 0 0 LopN2 0 0 LopN3 0 0;
    0 LopN1 0 0 LopN2 0 0 LopN3 0;
    0 0 0 0 0 0 0 0 0 ];
b = physics.bodyForce.magnitude*physics.bodyForce.direction;


duN = [ u(1)* dDotN(1)+u(4)*dDotN(4)+u(7)*dDotN(7), u(2)*dDotN(1)+u(5)*dDotN(4)+u(8)*dDotN(7), 0;
    u(1)* dDotN(2)+u(4)*dDotN(5)+u(7)*dDotN(8), u(2)*dDotN(2)+u(5)*dDotN(5)+u(8)*dDotN(8), 0;
    0,                                          0,                                         0 ];
CTangentBOperator = duN*N;

%% 4. Compute the diffusion matrix at the Gauss Point
De = physics.nue*(dN'*dN);

dDedN = [2 * dN1dx * u1x + dN2dx * u2x + dN3dx * u3x 2 * dN1dy * u1x + dN2dy * u2x + dN3dy * u3x dN1dx * u2x dN1dy * u2x dN1dx * u3x dN1dy * u3x; 
    2 * dN1dx * u1y + dN2dx * u2y + dN3dx * u3y 2 * dN1dy * u1y + dN2dy * u2y + dN3dy * u3y dN1dx * u2y dN1dy * u2y dN1dx * u3y dN1dy * u3y;
    0 0 0 0 0 0; 
    dN2dx * u1x dN2dy * u1x dN1dx * u1x + 2 * dN2dx * u2x + dN3dx * u3x dN1dy * u1x + 2 * dN2dy * u2x + dN3dy * u3x dN2dx * u3x dN2dy * u3x;
    dN2dx * u1y dN2dy * u1y dN1dx * u1y + 2 * dN2dx * u2y + dN3dx * u3y dN1dy * u1y + 2 * dN2dy * u2y + dN3dy * u3y dN2dx * u3y dN2dy * u3y; 
    0 0 0 0 0 0; 
    dN3dx * u1x dN3dy * u1x dN3dx * u2x dN3dy * u2x dN1dx * u1x + dN2dx * u2x + 2 * dN3dx * u3x dN1dy * u1x + dN2dy * u2x + 2 * dN3dy * u3x; 
    dN3dx * u1y dN3dy * u1y dN3dx * u2y dN3dy * u2y dN1dx * u1y + dN2dx * u2y + 2 * dN3dx * u3y dN1dy * u1y + dN2dy * u2y + 2 * dN3dy * u3y; 
    0 0 0 0 0 0;];
dDedN =  physics.nue*dDedN;
% Verification with finite difference


dDedX = dDedN*dNdX;

%% 5. Compute the convection matrix at the Gauss Point
Ce = N'*BoperatorConvection;

dCedN = [Gp1 * ux * u1x Gp1 * uy * u1x Gp1 * ux * u2x Gp1 * uy * u2x Gp1 * ux * u3x Gp1 * uy * u3x; 
    Gp1 * ux * u1y Gp1 * uy * u1y Gp1 * ux * u2y Gp1 * uy * u2y Gp1 * ux * u3y Gp1 * uy * u3y;
    0 0 0 0 0 0; 
    Gp2 * ux * u1x Gp2 * uy * u1x Gp2 * ux * u2x Gp2 * uy * u2x Gp2 * ux * u3x Gp2 * uy * u3x;
    Gp2 * ux * u1y Gp2 * uy * u1y Gp2 * ux * u2y Gp2 * uy * u2y Gp2 * ux * u3y Gp2 * uy * u3y; 
    0 0 0 0 0 0;
    Gp3 * ux * u1x Gp3 * uy * u1x Gp3 * ux * u2x Gp3 * uy * u2x Gp3 * ux * u3x Gp3 * uy * u3x;
    Gp3 * ux * u1y Gp3 * uy * u1y Gp3 * ux * u2y Gp3 * uy * u2y Gp3 * ux * u3y Gp3 * uy * u3y; 
    0 0 0 0 0 0;];
dCedX = dCedN * dNdX;


%% 6. Compute the convection-erro = analy(:,1) - FDtangent matrix at the Gauss Point
KTe = (N'*CTangentBOperator);
dKTedN = [u1x * Gp1 ^ 2 0 u2x * Gp1 ^ 2 0 u3x * Gp1 ^ 2 0;
          0 u1x * Gp1 ^ 2 0 u2x * Gp1 ^ 2 0 u3x * Gp1 ^ 2;
          0 0 0 0 0 0;
          Gp1 * u1x * Gp2 0 Gp1 * u2x * Gp2 0 Gp1 * u3x * Gp2 0;
          0 Gp1 * u1x * Gp2 0 Gp1 * u2x * Gp2 0 Gp1 * u3x * Gp2;
          0 0 0 0 0 0; 
          Gp1 * u1x * Gp3 0 Gp1 * u2x * Gp3 0 Gp1 * u3x * Gp3 0;
          0 Gp1 * u1x * Gp3 0 Gp1 * u2x * Gp3 0 Gp1 * u3x * Gp3;
          0 0 0 0 0 0;];
dKTedX = dKTedN * dNdX;


%% 7. Compute the velocity-pressure coupling matrix at the Gauss Point
Pe = -1*dDotN' * Np;

dPedN = [-0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3 0 0 0 0 0; 
          0 -0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3 0 0 0 0;
          0 0 0 0 0 0;
          0 0 -0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3 0 0 0;
          0 0 0 -0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3 0 0;
          0 0 0 0 0 0;
          0 0 0 0 -0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3 0;
          0 0 0 0 0 -0.1e1 * Gp1 * p1 - 0.1e1 * Gp2 * p2 - 0.1e1 * Gp3 * p3;
          0 0 0 0 0 0;];
dPedX = dPedN * dNdX;



%% Computing the variation of the pressure-velocity(Pe') matrix at Gauss Point
dPvdN = [0 0 0 0 0 0;
         0 0 0 0 0 0;
         -0.1e1 * Gp1 * u1x -0.1e1 * Gp1 * u1y -0.1e1 * Gp1 * u2x -0.1e1 * Gp1 * u2y -0.1e1 * Gp1 * u3x -0.1e1 * Gp1 * u3y; 
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         -0.1e1 * Gp2 * u1x -0.1e1 * Gp2 * u1y -0.1e1 * Gp2 * u2x -0.1e1 * Gp2 * u2y -0.1e1 * Gp2 * u3x -0.1e1 * Gp2 * u3y; 
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         -0.1e1 * Gp3 * u1x -0.1e1 * Gp3 * u1y -0.1e1 * Gp3 * u2x -0.1e1 * Gp3 * u2y -0.1e1 * Gp3 * u3x -0.1e1 * Gp3 * u3y;];
dPvdX = dPvdN * dNdX;


%% 8. Compute the stabilization matrices on the Gauss Point
%% 8i. Matrices related to diffusion nue*Laplacian(v)

% Diffusion-diffusion stabilization:
Dd = -1*(physics.nue^2)*(laplaceN'*laplaceN);
dDddN = zeros(9,6);
dDddX = dDddN * dNdX;

% Diffusion-convection stabilization:
Dc = physics.nue*(laplaceN'*BoperatorConvection);
dDcdN = zeros(9,6);
dDcdX = dDcdN *dNdX;

% Diffusion-pressure stabilization:
Dp = physics.nue*(laplaceN'*dNp);
dDpdN = zeros(9,6);
dDpdX = dDcdN *dNdX;

%% 8ii. Matrices related to convection v.Gradient(v)
% Convection-diffusion stabilization:
Cd = - Dc';
dCddN = zeros(9,6);
dCddX = dCddN * dNdX;


% Convection-convection stabilization:
Cc = BoperatorConvection'*BoperatorConvection;
dCcdN = [2 * (dN1dx * ux + dN1dy * uy) * u1x * ux + ux * (dN2dx * ux + dN2dy * uy) * u2x + ux * (dN3dx * ux + dN3dy * uy) * u3x 2 * (dN1dx * ux + dN1dy * uy) * u1x * uy + uy * (dN2dx * ux + dN2dy * uy) * u2x + uy * (dN3dx * ux + dN3dy * uy) * u3x (dN1dx * ux + dN1dy * uy) * ux * u2x (dN1dx * ux + dN1dy * uy) * uy * u2x (dN1dx * ux + dN1dy * uy) * ux * u3x (dN1dx * ux + dN1dy * uy) * uy * u3x;
         2 * (dN1dx * ux + dN1dy * uy) * u1y * ux + ux * (dN2dx * ux + dN2dy * uy) * u2y + ux * (dN3dx * ux + dN3dy * uy) * u3y 2 * (dN1dx * ux + dN1dy * uy) * u1y * uy + uy * (dN2dx * ux + dN2dy * uy) * u2y + uy * (dN3dx * ux + dN3dy * uy) * u3y (dN1dx * ux + dN1dy * uy) * ux * u2y (dN1dx * ux + dN1dy * uy) * uy * u2y (dN1dx * ux + dN1dy * uy) * ux * u3y (dN1dx * ux + dN1dy * uy) * uy * u3y;
         0 0 0 0 0 0;
         ux * (dN2dx * ux + dN2dy * uy) * u1x uy * (dN2dx * ux + dN2dy * uy) * u1x (dN1dx * ux + dN1dy * uy) * u1x * ux + 2 * ux * (dN2dx * ux + dN2dy * uy) * u2x + ux * (dN3dx * ux + dN3dy * uy) * u3x (dN1dx * ux + dN1dy * uy) * u1x * uy + 2 * uy * (dN2dx * ux + dN2dy * uy) * u2x + uy * (dN3dx * ux + dN3dy * uy) * u3x (dN2dx * ux + dN2dy * uy) * ux * u3x (dN2dx * ux + dN2dy * uy) * uy * u3x;
         ux * (dN2dx * ux + dN2dy * uy) * u1y uy * (dN2dx * ux + dN2dy * uy) * u1y (dN1dx * ux + dN1dy * uy) * u1y * ux + 2 * ux * (dN2dx * ux + dN2dy * uy) * u2y + ux * (dN3dx * ux + dN3dy * uy) * u3y (dN1dx * ux + dN1dy * uy) * u1y * uy + 2 * uy * (dN2dx * ux + dN2dy * uy) * u2y + uy * (dN3dx * ux + dN3dy * uy) * u3y (dN2dx * ux + dN2dy * uy) * ux * u3y (dN2dx * ux + dN2dy * uy) * uy * u3y;
         0 0 0 0 0 0;
         ux * (dN3dx * ux + dN3dy * uy) * u1x uy * (dN3dx * ux + dN3dy * uy) * u1x ux * (dN3dx * ux + dN3dy * uy) * u2x uy * (dN3dx * ux + dN3dy * uy) * u2x (dN1dx * ux + dN1dy * uy) * u1x * ux + ux * (dN2dx * ux + dN2dy * uy) * u2x + 2 * ux * (dN3dx * ux + dN3dy * uy) * u3x (dN1dx * ux + dN1dy * uy) * u1x * uy + uy * (dN2dx * ux + dN2dy * uy) * u2x + 2 * uy * (dN3dx * ux + dN3dy * uy) * u3x; 
         ux * (dN3dx * ux + dN3dy * uy) * u1y uy * (dN3dx * ux + dN3dy * uy) * u1y ux * (dN3dx * ux + dN3dy * uy) * u2y uy * (dN3dx * ux + dN3dy * uy) * u2y (dN1dx * ux + dN1dy * uy) * u1y * ux + ux * (dN2dx * ux + dN2dy * uy) * u2y + 2 * ux * (dN3dx * ux + dN3dy * uy) * u3y (dN1dx * ux + dN1dy * uy) * u1y * uy + uy * (dN2dx * ux + dN2dy * uy) * u2y + 2 * uy * (dN3dx * ux + dN3dy * uy) * u3y; 
         0 0 0 0 0 0;];
dCcdX = dCcdN * dNdX;


% Convection-pressure stabilization:
Cp = BoperatorConvection'*dNp;
dCpdN = [ux * dN1dx * p1 + (dN1dx * ux + dN1dy * uy) * p1 + ux * dN2dx * p2 + ux * dN3dx * p3 uy * dN1dx * p1 + uy * dN2dx * p2 + uy * dN3dx * p3 (dN1dx * ux + dN1dy * uy) * p2 0 (dN1dx * ux + dN1dy * uy) * p3 0;
         ux * dN1dy * p1 + ux * dN2dy * p2 + ux * dN3dy * p3 uy * dN1dy * p1 + (dN1dx * ux + dN1dy * uy) * p1 + uy * dN2dy * p2 + uy * dN3dy * p3 0 (dN1dx * ux + dN1dy * uy) * p2 0 (dN1dx * ux + dN1dy * uy) * p3;
         0 0 0 0 0 0;
         (dN2dx * ux + dN2dy * uy) * p1 0 ux * dN1dx * p1 + ux * dN2dx * p2 + (dN2dx * ux + dN2dy * uy) * p2 + ux * dN3dx * p3 uy * dN1dx * p1 + uy * dN2dx * p2 + uy * dN3dx * p3 (dN2dx * ux + dN2dy * uy) * p3 0;
         0 (dN2dx * ux + dN2dy * uy) * p1 ux * dN1dy * p1 + ux * dN2dy * p2 + ux * dN3dy * p3 uy * dN1dy * p1 + uy * dN2dy * p2 + (dN2dx * ux + dN2dy * uy) * p2 + uy * dN3dy * p3 0 (dN2dx * ux + dN2dy * uy) * p3;
         0 0 0 0 0 0; 
         (dN3dx * ux + dN3dy * uy) * p1 0 (dN3dx * ux + dN3dy * uy) * p2 0 ux * dN1dx * p1 + ux * dN2dx * p2 + ux * dN3dx * p3 + (dN3dx * ux + dN3dy * uy) * p3 uy * dN1dx * p1 + uy * dN2dx * p2 + uy * dN3dx * p3; 
         0 (dN3dx * ux + dN3dy * uy) * p1 0 (dN3dx * ux + dN3dy * uy) * p2 ux * dN1dy * p1 + ux * dN2dy * p2 + ux * dN3dy * p3 uy * dN1dy * p1 + uy * dN2dy * p2 + uy * dN3dy * p3 + (dN3dx * ux + dN3dy * uy) * p3;
         0 0 0 0 0 0;];
dCpdX = dCpdN * dNdX;

%% 8iii. Matrices related to pressure gradient Nabla(p)

% Pressure-diffusion stabilization: (Continuity equation)
Pd = - Dp';
dPddN = zeros(9,6);
dPddX =  dPddN * dNdX;


% Pressure-convection stabilization: (Continuity equation)
Pc = Cp';
dPcdN = [0 0 0 0 0 0; 0 0 0 0 0 0;
        ux * dN1dx * u1x + (dN1dx * ux + dN1dy * uy) * u1x + ux * dN1dy * u1y + (dN2dx * ux + dN2dy * uy) * u2x + (dN3dx * ux + dN3dy * uy) * u3x uy * dN1dx * u1x + uy * dN1dy * u1y + (dN1dx * ux + dN1dy * uy) * u1y + (dN2dx * ux + dN2dy * uy) * u2y + (dN3dx * ux + dN3dy * uy) * u3y ux * dN1dx * u2x + ux * dN1dy * u2y uy * dN1dx * u2x + uy * dN1dy * u2y ux * dN1dx * u3x + ux * dN1dy * u3y uy * dN1dx * u3x + uy * dN1dy * u3y;
        0 0 0 0 0 0; 
        0 0 0 0 0 0; 
        ux * dN2dx * u1x + ux * dN2dy * u1y uy * dN2dx * u1x + uy * dN2dy * u1y (dN1dx * ux + dN1dy * uy) * u1x + ux * dN2dx * u2x + (dN2dx * ux + dN2dy * uy) * u2x + ux * dN2dy * u2y + (dN3dx * ux + dN3dy * uy) * u3x (dN1dx * ux + dN1dy * uy) * u1y + uy * dN2dx * u2x + uy * dN2dy * u2y + (dN2dx * ux + dN2dy * uy) * u2y + (dN3dx * ux + dN3dy * uy) * u3y ux * dN2dx * u3x + ux * dN2dy * u3y uy * dN2dx * u3x + uy * dN2dy * u3y;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        ux * dN3dx * u1x + ux * dN3dy * u1y uy * dN3dx * u1x + uy * dN3dy * u1y ux * dN3dx * u2x + ux * dN3dy * u2y uy * dN3dx * u2x + uy * dN3dy * u2y (dN1dx * ux + dN1dy * uy) * u1x + (dN2dx * ux + dN2dy * uy) * u2x + ux * dN3dx * u3x + (dN3dx * ux + dN3dy * uy) * u3x + ux * dN3dy * u3y (dN1dx * ux + dN1dy * uy) * u1y + (dN2dx * ux + dN2dy * uy) * u2y + uy * dN3dx * u3x + uy * dN3dy * u3y + (dN3dx * ux + dN3dy * uy) * u3y;];

dPcdX = dPcdN * dNdX;

% Pressure-pressure stabilization: (Continuity equation)
Pp = dNp'*dNp;
dPpdN = [0 0 0 0 0 0; 0 0 0 0 0 0; 
        2 * dN1dx * p1 + p2 * dN2dx + p3 * dN3dx 2 * dN1dy * p1 + dN2dy * p2 + dN3dy * p3 dN1dx * p2 dN1dy * p2 dN1dx * p3 dN1dy * p3;
        0 0 0 0 0 0; 0 0 0 0 0 0;
        dN2dx * p1 dN2dy * p1 dN1dx * p1 + 2 * p2 * dN2dx + p3 * dN3dx dN1dy * p1 + 2 * dN2dy * p2 + dN3dy * p3 dN2dx * p3 dN2dy * p3;
        0 0 0 0 0 0; 0 0 0 0 0 0;
        dN3dx * p1 dN3dy * p1 dN3dx * p2 dN3dy * p2 dN1dx * p1 + p2 * dN2dx + 2 * p3 * dN3dx dN1dy * p1 + dN2dy * p2 + 2 * dN3dy * p3;];
dPpdX = dPpdN * dNdX;


%% 8iv. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
Div = dDotN'*dDotN; % (Continuity equation)
dDivdN = [2 * dN1dx * u1x + dN1dy * u1y + dN2dx * u2x + dN2dy * u2y + dN3dx * u3x + dN3dy * u3y dN1dx * u1y dN1dx * u2x dN1dx * u2y dN1dx * u3x dN1dx * u3y; 
          dN1dy * u1x dN1dx * u1x + 2 * dN1dy * u1y + dN2dx * u2x + dN2dy * u2y + dN3dx * u3x + dN3dy * u3y dN1dy * u2x dN1dy * u2y dN1dy * u3x dN1dy * u3y; 
          0 0 0 0 0 0;
          dN2dx * u1x dN2dx * u1y dN1dx * u1x + dN1dy * u1y + 2 * dN2dx * u2x + dN2dy * u2y + dN3dx * u3x + dN3dy * u3y dN2dx * u2y dN2dx * u3x dN2dx * u3y; 
          dN2dy * u1x dN2dy * u1y dN2dy * u2x dN1dx * u1x + dN1dy * u1y + dN2dx * u2x + 2 * dN2dy * u2y + dN3dx * u3x + dN3dy * u3y dN2dy * u3x dN2dy * u3y;
          0 0 0 0 0 0; 
          dN3dx * u1x dN3dx * u1y dN3dx * u2x dN3dx * u2y dN1dx * u1x + dN1dy * u1y + dN2dx * u2x + dN2dy * u2y + 2 * dN3dx * u3x + dN3dy * u3y dN3dx * u3y; 
          dN3dy * u1x dN3dy * u1y dN3dy * u2x dN3dy * u2y dN3dy * u3x dN1dx * u1x + dN1dy * u1y + dN2dx * u2x + dN2dy * u2y + dN3dx * u3x + 2 * dN3dy * u3y;
          0 0 0 0 0 0;];
dDivdX = dDivdN * dNdX;


%% 9. Compute the mass matrix
M = N'*N;
dMdN = zeros(9,6); 
dMdX = dMdN * dNdX;

%% 10. Compute the stabilization matrices on the mass
%% 10i. Mass-diffusion stabilization
Md = physics.nue*laplaceN'*N;
dMddN = zeros(9,6);
dMddX = dMdN * dNdX;

%% 10ii. Mass-convection stabilization
Mc = BoperatorConvection'*N;
dMcdN = [Gp1 * ux * u1x + Gp2 * ux * u2x + Gp3 * ux * u3x Gp1 * uy * u1x + Gp2 * uy * u2x + Gp3 * uy * u3x 0 0 0 0; 
         Gp1 * ux * u1y + Gp2 * ux * u2y + Gp3 * ux * u3y Gp1 * uy * u1y + Gp2 * uy * u2y + Gp3 * uy * u3y 0 0 0 0;
         0 0 0 0 0 0;
         0 0 Gp1 * ux * u1x + Gp2 * ux * u2x + Gp3 * ux * u3x Gp1 * uy * u1x + Gp2 * uy * u2x + Gp3 * uy * u3x 0 0;
         0 0 Gp1 * ux * u1y + Gp2 * ux * u2y + Gp3 * ux * u3y Gp1 * uy * u1y + Gp2 * uy * u2y + Gp3 * uy * u3y 0 0;
         0 0 0 0 0 0;
         0 0 0 0 Gp1 * ux * u1x + Gp2 * ux * u2x + Gp3 * ux * u3x Gp1 * uy * u1x + Gp2 * uy * u2x + Gp3 * uy * u3x;
         0 0 0 0 Gp1 * ux * u1y + Gp2 * ux * u2y + Gp3 * ux * u3y Gp1 * uy * u1y + Gp2 * uy * u2y + Gp3 * uy * u3y;
         0 0 0 0 0 0;];
dMcdX = dMcdN * dNdX;

%% 10iii. Mass-pressure stabilization
Mp = dNp'*N;
dMpdN = [0 0 0 0 0 0; 0 0 0 0 0 0;
         Gp1 * u1x + Gp2 * u2x + Gp3 * u3x Gp1 * u1y + Gp2 * u2y + Gp3 * u3y 0 0 0 0; 
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 Gp1 * u1x + Gp2 * u2x + Gp3 * u3x Gp1 * u1y + Gp2 * u2y + Gp3 * u3y 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 Gp1 * u1x + Gp2 * u2x + Gp3 * u3x Gp1 * u1y + Gp2 * u2y + Gp3 * u3y;];
dMpdX = dMpdN * dNdX;


%% 11. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
% Stabilization vector related to nue*Laplacian(v) (reaction):
Fd = physics.nue*laplaceN';
dFdN = zeros(9,6); 
dFdX = dFdN * dNdX;

% Stabilization due to pressure gradient Nabla(p): (Continuity equation)
Fp = dNp';
dFpdN = zeros(9,6);
dFpdX = dFdN * dNdX;

% 11ii Convection-Body force stabilization
Fc = -1*(BoperatorConvection');
dFcdN = zeros(9,6); 
dFcdX = dFcdN * dNdX;

%% Anal till here

analy = dDedX + dCedX + dPedX  + tauM*(dDddX + dDcdX + dDpdX + dCddX + dCcdX + dCpdX + dPddX + dPcdX + dPpdX) ...
                                      + ( (Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp)*u )* dtauMdX'  ...
                                      + tauC*dDivdX ...
                                      + (Div*u)*dtauCdX';
                                  
manaly = dMdX + tauM*(dMddX + dMcdX + dMpdX) ...
             + ((Md + Mc + Mp)*u)*dtauMdX';                            
% analy = tauC*dDivdN;
% analy =  ( (Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp)*u )* dtauMdX';
erro = 100*(analy(:,2) - FD)/max(abs(FD))
merro = 100*(manaly(:,2) - FDM)/max(abs(FDM))




%% 12. Compute the element stiffness matrix
% Linear stiffness matrix needed for the residual computation
Ke = De + Ce + Pe - Pe' + tauM*(Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp) + tauC*Div;
dKedX = dDedX + dCedX + dPedX - dPvdX + tauM*(dDddX + dDcdX + dDpdX + dCddX + dCcdX + dCpdX + dPddX + dPcdX + dPpdX) ...
                                      + ( (Dd + Dc + Dp + Cd + Cc + Cp + Pd + Pc + Pp)*u )* dtauMdX'  ...
                                      + tauC*dDivdX ...
                                      + (Div*u)*dtauCdX';

%% 13. Compute the element stabilized mass matrix
Me = M + tauM*(Md + Mc + Mp);
dMedX = dMdX + tauM*(dMddX + dMcdX + dMpdX) ...
             + ((Md + Mc + Mp)*u)*dtauMdX';

%% 14. Compute the stabilized force vector due to source at the Gauss point

% Compute the complete right-hand side vector at the Gauss Point
% Compute the force vector at the Gauss Point
Fe = (N' + tauM*(Fd + Fp + Fc))*b;
dFedX = zeros(9,6);

% End of function
end
