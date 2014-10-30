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
function B = computeBOperatorMatrixForPlateInMembraneActionCST(nodes)
%% Function documentation
%
% Returns the B-operator matrix at the element level for the case of
% discretizing a plate in membrane action problem using the Constant Strain
% Triangle (CST)
%
%   Input :
%   nodes : The nodes of the CST in a counterclock-wise fashion
% 
%  Output :
%       B : The B-operator matrix in the element level
%
% Function layout:
%
% 0. Read input
%
% 1. Compute the B-operator matrix in the element level
%
%% Function main body

%% 0. Read input

% Get the nodes of the mesh in a counter-clockwise fashion
%
%        i ------------- k
%          \           /
%           \         /
%            \       /
%             \     /  
%              \   /
%               \ /
%                j
vertix_i = nodes(1,:);
vertix_j = nodes(2,:);
vertix_k = nodes(3,:);

% Compute the permutation symbols
% For basis function Ni:
% yjk:
yjk = vertix_j(2)-vertix_k(2);
% xkj:
xkj = vertix_k(1) - vertix_j(1);

% For basis function Nj:
% yik:
yik = -(vertix_i(2)-vertix_k(2));
% xki:
xki = -(vertix_k(1) - vertix_i(1));

% For basis function Nj:
% yij:
yij = vertix_i(2)-vertix_j(2);
% xji:
xji = vertix_j(1) - vertix_i(1);

% The area of the triangle
% Delta = det(nodes)/2;
Delta = abs((nodes(1,1)-nodes(3,1))*(nodes(2,2)-nodes(1,2))-(nodes(1,1)-nodes(2,1))*(nodes(3,2)-nodes(1,2)))/2;

%% 1. Compute the B-operator matrix in the element level

preFactor_B = 1/2/Delta;
B = preFactor_B*[yjk 0 yik 0 yij 0
                 0 xkj 0 xki 0 xji
                 xkj yjk xki yik xji yij];

end

