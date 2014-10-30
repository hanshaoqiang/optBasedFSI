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
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = load_line_triangle_2d(p,U,q,V,CP,F_old,mesh,p,ub,vb)
%% Function documentation
%
% Returns the consistent load vector corresponding to a boundary load of
% magnitude p and a 2D mesh consisting of triangular elements
%
%   Input :
%     p,q : polynomial degrees of the NURBS geometry
%     U,V : knot vectors in u,v-direction
%      CP : set of control point coordinates and weights
%   F_old : Existing load vector
%    mesh : The triangular mesh
%       p : load magnitude [N/m] (line load->pressure)
%      ub : u-extension of the load area
%      vb : v-extension of the load area
%
% Function layout :
%
%
%% Function main body




end

