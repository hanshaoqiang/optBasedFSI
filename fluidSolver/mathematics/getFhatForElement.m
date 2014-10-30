function [ Fhat ] = getFhatForElement( nodes, dt, u)
%   This method contains routene for calculating the Gradient matrix for
%   the ALE mapping which is called here Fhat
%   
%   Input:
%       nodes   -- Matrix containing the coordinates of nodes of the current tringular
%                   element. Its is a 3x3 matrix.
%       dt      -- Time step of the simulation.
%       u       -- velocity vector containing the x and y velocities at the
%                  nodes.
%   Output:
%       Fhat    -- Derivative matrix for the ALE mapping.

% Obtaining the x and y coordinates of all the nodes.
% NODE 1
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

% Zeros matrix for making Fhat
z = zeros(3);

%% Fromulating the necessary deravatives for Fhat
% L stand for Lambda in the following formulation
% Terms of first node
u1x = u(1);
u1y = u(2);
dL1x_dx = 1 + dt*dN1_dx*u1x;
dL1x_dy = 1 + dt*dN1_dy*u1x;

dL1y_dx = 1 + dt*dN1_dx*u1y;
dL1y_dy = 1 + dt*dN1_dy*u1y;

m1 = [dL1x_dx dL1x_dy;
      dL1y_dx dL1y_dy ];
m1 = m1^-1;
M1 = [m1(1,1) m1(1,2) 0;
      m1(2,1) m1(2,2) 0;
      0       0       0 ];

% Terms for Second Node
u2x = u(3);
u2y = u(4);
dL2x_dx = 1 + dt*dN2_dx*u2x;
dL2x_dy = 1 + dt*dN2_dy*u2x;

dL2y_dx = 1 + dt*dN2_dx*u2y;
dL2y_dy = 1 + dt*dN2_dy*u2y;

m2 = [dL2x_dx dL2x_dy;
      dL2y_dx dL2y_dy ];
m2 = m2^-1;
M2 = [m2(1,1) m2(1,2) 0;
      m2(2,1) m2(2,2) 0;
      0       0       0 ];

% Terms for Third node
u3x = u(5);
u3y = u(6);

dL3x_dx = 1 + dt*dN3_dx*u3x;
dL3x_dy = 1 + dt*dN3_dy*u3x;

dL3y_dx = 1 + dt*dN3_dx*u3y;
dL3y_dy = 1 + dt*dN3_dy*u3y;

m3 = [dL3x_dx dL3x_dy;
      dL3y_dx dL3y_dy ];

m3 = m3^-1;
M3 = [m3(1,1) m3(1,2) 0;
      m3(2,1) m3(2,2) 0;
      0       0       0 ];


% Finally comming to Fhat
% It should be a 9x9 matrix
Fhat = [M1 z z;
        z M2 z;
        z z  M3];
    
    
% End of the function
end