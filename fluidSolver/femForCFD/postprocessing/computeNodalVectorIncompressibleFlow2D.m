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
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unknownVector = computeNodalVectorIncompressibleFlow2D(i,p,u,U,j,q,v,V,CP,elementNodalVector)
%% Function documentation
%
% Returns the element nodal vector vector u = [ux uy p]' at the given 
% parametric location for a 2D incompressible flow problem.
%
%         Input : 
%           p,q : polynomial degrees
%           U,V : knot vectors in u,v-direction
%            CP : Control Point coordinates and weights
%           u,v : Parametric coordinate of the NURBS surface
%        scalar : element transport vector
%
%        Output :
% unknownVector : The unknown vector u = [ux uy p]' at the given parametric
%                 location
%
% Function layout :
%
% 1. Compute the shape function matrix
%
% 2. Compute the unknown vector as a product u = Rmatrix*d
%
%% Function main body

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 3*ne;

% Compute NURBS basis functions
R = computeNurbsBasisFunctions2D(i,p,u,U,j,q,v,V,CP);

% Compute the shape function matrix
Rmatrix = zeros(3,ndof);

for i=1:ne
    Rmatrix(1,3*i-2) = R(i);
    Rmatrix(2,3*i-1) = R(i);
    Rmatrix(3,3*i) = R(i);
end

% displacement
unknownVector = Rmatrix*elementNodalVector;

end