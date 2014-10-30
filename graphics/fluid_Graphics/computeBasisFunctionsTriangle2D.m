%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Dr.-Ing. Roland Wüchner                                               %
%   Dipl.-Math. Andreas Apostolatos (andreas.apostolatos@tum.de)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = computeBasisFunctionsTriangle2D(vertix_i,vertix_j,vertix_k,x,y)
%% Function documentation
%
% Returns the triangular shape functions given the three vertices of the
% triangle and the location on the x-y plane. The vertices must be provided
% in a counterclock-wise fashion:
%
%               k
%              / \
%             /   \
%            /     \ 
%           /       \
%          /         \
%         /           \
%        i-------------j
%
%     Input :
%    vertix : The coordinates of the vertices counterclockwise
%       x,y : Location on the x,y plane to compute the shape functions
%
%    Output :
%         N : The evaluated basis functions at x,y: N=[Ni Nj Nk]
%
% Function layout :
%
% 1. Compute the area of the triangle
%
% 2. Compute the permutations
%
% 3. Compute the basis functions for the linear triangle at (x,y)
%
%% Function main body

%% 1. Compute the area of the triangle

Delta = .5*det([1 vertix_i(1) vertix_i(2);
             1 vertix_j(1) vertix_j(2);
             1 vertix_k(1) vertix_k(2)]);

%% 2. Compute the permutations

% For basis function Ni:
% zi:
zi = vertix_j(1)*vertix_k(2)-vertix_k(1)*vertix_j(2);
% yjk:
yjk = vertix_j(2)-vertix_k(2);
% xkj:
xkj = vertix_k(1) - vertix_j(1);

% For basis function Nj:
% zj:
zj = (vertix_k(1)*vertix_i(2)-vertix_i(1)*vertix_k(2));
% yik:
yik = -(vertix_i(2)-vertix_k(2));
% xki:
xki = -(vertix_k(1) - vertix_i(1));

% For basis function Nj:
% zk:
zk = vertix_i(1)*vertix_j(2)-vertix_j(1)*vertix_i(2);
% yij:
yij = vertix_i(2)-vertix_j(2);
% xji:
xji = vertix_j(1) - vertix_i(1);

%% 3. Compute the basis functions for the linear triangle at (x,y)

% Ni:
Ni = (zi+yjk*x+xkj*y)/2/Delta;

% Nj:
Nj = (zj+yik*x+xki*y)/2/Delta;

% Nk:
Nk = (zk+yij*x+xji*y)/2/Delta;

% vector containing all the basis functions
N = [Ni Nj Nk];

end

