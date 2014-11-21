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
%   Reza Najian (M.Sc)                 (reza.najian-asl@tum.de)           %
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de          %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dNdx ] = getVariationOfVelocityShapeFunctionsWithCord( GP, nodes )
%   This method contains routene for calculating the derivatives of Nodal
%   basis on a particular element.
%   
%   Input:
%       nodes   -- Matrix containing the coordinates of nodes of the current tringular
%                   element. Its is a 3x3 matrix.
%       GP      -- Gauss point where the matrices are needed.
%   Output:
%       dNdx    -- The variation of the Shape functions with respect to the
%                  Coordinates (6x6 matrix) 


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
area_unmod = ( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );

%% Calculating all the derivatives of the Shape functions
% NOTE : Basis Function used are of Order 1
% #############

% variation of dN1dx
dN1dx_x1 = (-1/8/(area^3))*((y2-y3)^2)*area_unmod;
% dN1dx_x1 = (y2 - y3)^2 * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN1dx_x2 = (-1/8/(area^3))*(y2-y3)*(y3-y1)*area_unmod;
dN1dx_x3 = (-1/8/(area^3))*((y2-y3)*(y1-y2))*area_unmod;
dN1dx_y1 = (-1/8/(area^3))*((y2-y3)*(x3-x2))*area_unmod;
dN1dx_y2 = (1/2/area) + (-1/8/(area^3))*((y2-y3)*(x1-x3))*area_unmod;
dN1dx_y3 = (-1/2/area) + (-1/8/(area^3))*((y2-y3)*(x2-x1))*area_unmod;
% 
% % variation of dN1dy
dN1dy_x1 = (-1/8/(area^3))*((x3-x2)*(y2-y3))*area_unmod;
% dN1dy_x1 = (y2 - y3)*(x3-x2) * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN1dy_x2 = (-1/2/area) + (-1/8/(area^3))*((x3-x2)*(y3-y1))*area_unmod;
dN1dy_x3 = (1/2/area) + (-1/8/(area^3))*((x3-x2)*(y1-y2))*area_unmod;
dN1dy_y1 = (-1/8/(area^3))*((x3-x2)*(x3-x2))*area_unmod;
dN1dy_y2 = (-1/8/(area^3))*((x3-x2)*(x1-x3))*area_unmod;
dN1dy_y3 = (-1/8/(area^3))*((x3-x2)*(x2-x1))*area_unmod;
% 
% % variation of dN2dx
dN2dx_x1 = (-1/8/(area^3))*((y3-y1)*(y2-y3))*area_unmod;
% dN2dx_x1 = (y3 - y1)*(y2-y3) * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN2dx_x2 = (-1/8/(area^3))*((y3-y1)*(y3-y1))*area_unmod;
dN2dx_x3 = (-1/8/(area^3))*((y3-y1)*(y1-y2))*area_unmod;
dN2dx_y1 = (-1/2/area) + (-1/8/(area^3))*((y3-y1)*(x3-x2))*area_unmod;
dN2dx_y2 = (-1/8/(area^3))*((y3-y1)*(x1-x3))*area_unmod;
dN2dx_y3 = (1/2/area) + (-1/8/(area^3))*((y3-y1)*(x2-x1))*area_unmod;
% 
% % variation of dN2dy
dN2dy_x1 = (1/2/area) + (-1/8/(area^3))*((x1-x3)*(y2-y3))*area_unmod;
% dN2dy_x1 = 1/2/area + (x1-x3)*(y2-y3) * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN2dy_x2 = (-1/8/(area^3))*((x1-x3)*(y3-y1))*area_unmod;
dN2dy_x3 = (-1/2/area) + (-1/8/(area^3))*((x1-x3)*(y1-y2))*area_unmod;
dN2dy_y1 = (-1/8/(area^3))*((x1-x3)*(x3-x2))*area_unmod;
dN2dy_y2 = (-1/8/(area^3))*((x1-x3)*(x1-x3))*area_unmod;
dN2dy_y3 = (-1/8/(area^3))*((x1-x3)*(x2-x1))*area_unmod;

% % variation of dN3dx
dN3dx_x1 = (-1/8/(area^3))*((y1-y2)*(y2-y3))*area_unmod;
% dN3dx_x1 = (y1-y2)*(y2-y3) * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN3dx_x2 = (-1/8/(area^3))*((y1-y2)*(y3-y1))*area_unmod;
dN3dx_x3 = (-1/8/(area^3))*((y1-y2)*(y1-y2))*area_unmod;
dN3dx_y1 = (1/2/area) + (-1/8/(area^3))*((y1-y2)*(x3-x2))*area_unmod;
dN3dx_y2 = (-1/2/area) + (-1/8/(area^3))*((y1-y2)*(x1-x3))*area_unmod;
dN3dx_y3 = (-1/8/(area^3))*((y1-y2)*(x2-x1))*area_unmod;
% 
% % variation of dN3dy
% dN3dy_x1 = (-1/2/area) + (-1/8/(area^3))*((x2-x1)*(y2-y3))*area_unmod;
dN3dy_x1 = -1/2/area + (x2-x1)*(y2-y3) * (-1.0/4.0) * (1/area^2) * sign(area_unmod);
dN3dy_x2 = (1/2/area) + (-1/8/(area^3))*((x2-x1)*(y3-y1))*area_unmod;
dN3dy_x3 = (-1/8/(area^3))*((x2-x1)*(y1-y2))*area_unmod;
dN3dy_y1 = (-1/8/(area^3))*((x2-x1)*(x3-x2))*area_unmod;
dN3dy_y2 = (-1/8/(area^3))*((x2-x1)*(x1-x3))*area_unmod;
dN3dy_y3 = (-1/8/(area^3))*((x2-x1)*(x2-x1))*area_unmod;
% 
% 
% 
% % Constructing the variation matrix (6x6)
dNdx = [dN1dx_x1 dN1dx_y1 dN1dx_x2 dN1dx_y2 dN1dx_x3 dN1dx_y3;
        dN1dy_x1 dN1dy_y1 dN1dy_x2 dN1dy_y2 dN1dy_x3 dN1dy_y3;
        dN2dx_x1 dN2dx_y1 dN2dx_x2 dN2dx_y2 dN2dx_x3 dN2dx_y3;
        dN2dy_x1 dN2dy_y1 dN2dy_x2 dN2dy_y2 dN2dy_x3 dN2dy_y3;
        dN3dx_x1 dN3dx_y1 dN3dx_x2 dN3dx_y2 dN3dx_x3 dN3dx_y3;
        dN3dy_x1 dN3dy_y1 dN3dy_x2 dN3dy_y2 dN3dy_x3 dN3dy_y3;];
    

%% Verification
val_b = 0;
val = 0;
% Differential of Basis Function 1
dN1_dx_b = (y2-y3)/(2*area);
val_b(1) = dN1_dx_b;
dN1_dy_b = (x3-x2)/(2*area);
val_b(2) = dN1_dy_b;
% % Differential of Basis Function 2
dN2_dx_b = (y3-y1)/(2*area);
val_b(3) = dN2_dx_b;
dN2_dy_b = (x1-x3)/(2*area);
val_b(4) = dN2_dy_b;
% 
% Differential of Basis Function 3
dN3_dx_b = (y1-y2)/(2*area);
val_b(5) = dN3_dx_b;
dN3_dy_b = (x2-x1)/(2*area);
val_b(6) = dN3_dy_b;



delta = 1E-9;
% Peturbing x1
x2 = x2 + delta;
area = 0.5*abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
    
% Differential of Basis Function 1
dN1_dx = (y2-y3)/(2*area);
val(1) = dN1_dx;
dN1_dy = (x3-x2)/(2*area);
val(2) = dN1_dy;
% % Differential of Basis Function 2
dN2_dx = (y3-y1)/(2*area);
val(3) = dN2_dx;
dN2_dy = (x1-x3)/(2*area);
val(4) = dN2_dy;
% 
% Differential of Basis Function 3
dN3_dx = (y1-y2)/(2*area);
val(5) = dN3_dx;
dN3_dy = (x2-x1)/(2*area);
val(6) = dN3_dy;

% variation
v = (val - val_b)/delta;
x2 = x2- delta;
v = v';
v_a = dNdx(:,3);
error = v - v_a
% End of the Function        
end

