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
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates 
%       df
%      -----
%       dz
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
function [ df_dz ] = obtainVariationOfObjWithGeometry(nodes, objectiveFunction)

% Nodal coordinates
x1 = nodes(1,1); x2 = nodes(2,1);
y1 = nodes(1,2); y2 = nodes(2,2);

% Length of the given element
length = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

dx = x1 - x2;
dy = y1 - y2;

% Values of the objective function at the given element nodes.
f1 = objectiveFunction(1);
f2 = objectiveFunction(2);

% Calculating the the variation.
if(dx ~= 0)
    df_dx = f1/dx - f2/dx;
    df_dx = f1/length - f2/length;
else
    df_dx = 0;
    df_dx = f1/length - f2/length;
end
if(dy ~= 0)
    df_dy = f1/dy - f2/dy;
    df_dy = f1/length - f2/length;
else
    df_dy = 0;
    df_dy = f1/length - f2/length;
end

% Formulating the variance vector
df_dz =        [df_dx;
                df_dy;
                df_dx;
                df_dy];


% End of the function
end

