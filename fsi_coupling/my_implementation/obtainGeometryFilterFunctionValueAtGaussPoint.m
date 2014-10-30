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
% This method calculates and sends the value of the filter functions at the
% given gauss point. Since this is a 2D application, there will be edges on
% the FSI interface, so there filter functions are linear hat functions
% over the element. These can be choosen to be any thing. 
% 
%   Input:
%       GP          : Gaussian Coordinates for a Gauss point.
%       nodes       : Coordinates of the nodes forming the element.
%
%   Output:
%       Fcord       : This the matrix containing the values of the Filter
%                     function at the requested Gaussian point.
%   
function [ Fcord ] = obtainGeometryFilterFunctionValueAtGaussPoint(GP, nodes)
%%% Finding the values of the shape functions at gauss points
% Nodal coordinates
x1 = nodes(1,1); x2 = nodes(2,1);
y1 = nodes(1,2); y2 = nodes(2,2);
% Length of the given element
length = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

dx = x1 - x2;
dy = y1 - y2;
dx = length;
dy = length;

if(dx ~= 0)
    N1_x = 0.5 - (GP(1)/dx)*abs((x1+x2)/2 - x1);
    N2_x = 0.5 + (GP(1)/dx)*abs((x1+x2)/2 - x1);
else
    N1_x = 0;
    N2_x = 0;
    N1_x = 0.5 - (GP(1)/dx)*abs((x1+x2)/2 - x1);
    N2_x = 0.5 + (GP(1)/dx)*abs((x1+x2)/2 - x1);
end


if(dy ~= 0)
    N1_y = 0.5 - (GP(1)/dy)*abs((y1+y2)/2 - y1);
    N2_y = 0.5 + (GP(1)/dy)*abs((y1+y2)/2 - y1);
else
    N1_y = 0;
    N2_y = 0;
    N1_y = 0.5 - (GP(1)/dy)*abs((y1+y2)/2 - y1);
    N2_y = 0.5 + (GP(1)/dy)*abs((y1+y2)/2 - y1);
end

Fcord = [N1_x 0 N2_x 0;
         0 N1_y 0 N2_y];
     
% Fcord = [N1_x N1_y N2_x N2_y];

end

