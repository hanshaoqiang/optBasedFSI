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
function unknownVector = computeNodalVectorIncompressibleFlow2DGivenBasis(R,p,q,elementNodalVector)
%% Function documentation
%
% Returns the element nodal vector vector u = [ux uy p]' at the given 
% parametric location for a 2D incompressible flow problem.
%
%              Input :
%                  R : The pre-evaluated NURBS basis functions
%                p,q : polynomial degrees
% elementNodalVector : element transport vector
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