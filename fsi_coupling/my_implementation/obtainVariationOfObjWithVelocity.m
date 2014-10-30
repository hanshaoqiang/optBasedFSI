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
%   Aditya Ghantasala (M.Sc)           (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates 
%         df
%      --------
%       d(vel)
%
% This is used in the optimization based interface update.
% 
%   Input:
%       nodes               : Nodal coordinates making up the current element.
%       objectiveFunction   : The function field for which the variation is
%                             needed. This is provided only element wise.
%                             No need to provide it on the whole FSI
%                             interface.
%
%   Output:
%       df_dz               : The variation of the provided objective
%                             function field with respect to the geometry
%
function [ df_dvel ] = obtainVariationOfObjWithVelocity( velocity, nodes, objectiveFunction )
% Nodal coordinates
x1 = nodes(1,1); x2 = nodes(2,1);
y1 = nodes(1,2); y2 = nodes(2,2);

dx = x1 - x2;
dy = y1 - y2;

% Length of the given element
length = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
dx = length;
dy = length;


% velocities
vx1 = velocity(1,1); vx2 = velocity(2,1);
vy1 = velocity(1,2); vy2 = velocity(2,2);

% Variation of velocities
dv_dx = vx1/dx - vx2/dx;
dv_dy = vy1/dy - vy2/dy;

% Length of the given element
length_s = sqrt( (vx1-vx2)*(vx1-vx2) + (vy1-vy2)*(vy1-vy2) );
% Values of the objective function at the given element nodes.
f1 = objectiveFunction(1);
f2 = objectiveFunction(2);

dN1f_dx = 1/dx;
dN2f_dx = -1/dx;

dN1f_dy = 1/dy;
dN2f_dy = -1/dy;

% Calculating the the variation.
if( dx ~= 0)
    %df_dvx = f1/(vx1-vx2) - f2/(vx1-vx2);
    df_dvx = f1*dN1f_dx/dv_dx + f2*dN2f_dx/dv_dx;
else
    df_dvx = 0;
    df_dvx = f1*dN1f_dx/dv_dx + f2*dN2f_dx/dv_dx;
end

if( dy ~= 0)
    %df_dvy = f1/(vy1-vy2) - f2/(vy1-vy2);
    df_dvy = f1*dN1f_dy/dv_dy + f2*dN2f_dy/dv_dy;
else
    df_dvy = 0;
    df_dvy = f1*dN1f_dy/dv_dy + f2*dN2f_dy/dv_dy;
end

%df_dvel = [f1*(vx1-vx2)/length + f2*(vx2-vx1)/length 0;
%           0 f1*(vy1-vy2)/length + f2*(vy2-vy1)/length];
df_dvel = [ df_dvx ;
            df_dvy ;
            df_dvx ;
            df_dvy ];
        
        
if(norm(df_dvel) > 2)
%    fprintf('HELOOOOOO .. !!') ;
end

% End of the function
end

