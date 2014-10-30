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
function dN = computeCST2DBasisFunctionsAndFirstDerivatives(vertix_i,vertix_j,vertix_k,x,y)
%% Function documentation
%
% Returns the triangular shape functions and their derivatives with respect 
% to x and y coordinates given the three vertices of the triangle and the 
% location on the x-y plane. The vertices must be provided in a 
% counterclock-wise fashion:
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
%        dN : The evaluated basis functions and their first derivatives 
%             at x,y: dN=[Ni dNi/dx dNi/dy
%                         Nj dNj/dx dNj/dy
%                         Nk dNk/dx dNk/dy]
%
% Function layout :
%
% 1. Compute the area of the triangle
%
% 2. Compute the permutations
%
% 3. Compute the basis functions for the linear triangle at (x,y)
%
% 4. Compute the derivatives of the basis functions for the linear triangle w.r.t. to x at (x,y)
%
% 5. Compute the derivatives of the basis functions for the linear triangle w.r.t. to y at (x,y)
%
% 6. Assemble to the vector containing all the basis functions and their derivatives at (x,y)
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

%% 4. Compute the derivatives of the basis functions for the linear triangle w.r.t. to x at (x,y)

% Ni:
dNidx = yjk/2/Delta;

% Nj:
dNjdx = yik/2/Delta;

% Nk:
dNkdx = yij/2/Delta;

%% 5. Compute the derivatives of the basis functions for the linear triangle w.r.t. to y at (x,y)

% Ni:
dNidy = xkj/2/Delta;

% Nj:
dNjdy = xki/2/Delta;

% Nk:
dNkdy = xji/2/Delta;

%% 6. Assemble to the vector containing all the basis functions and their derivatives at (x,y)
dN=[Ni dNidx dNidy
    Nj dNjdx dNjdy
    Nk dNkdx dNkdy];

end

