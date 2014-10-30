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
function displacement = computeVectorTransport2DGivenBasis(R,p,q,transportVector)
%% Function documentation
%
% Returns the transport vector u = [ux uy]' at the given parametric
% location given the basis function
%
%        Input : 
%            R : The precomputed basis functions
%          p,q : polynomial degrees
%          u,v : Parametric coordinate of the NURBS surface
%       scalar : element transport vector
%
%       Output :
% displacement : The transport vector u = [ux uy]' at the given parametric
%                location
%
% Function layout :
%
% 1. Compute the shape function matrix
% 2. Compute the displacement vector as a product u = N*d
%
%% Function main body

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 2*ne;

% Compute the shape function matrix
N = zeros(2,ndof);

for i=1:ne
    N(1,2*i-1) = R(i);
    N(2,2*i) = R(i);
end

% displacement
displacement = N*transportVector;

end