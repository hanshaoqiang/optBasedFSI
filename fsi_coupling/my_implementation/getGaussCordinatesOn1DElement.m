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
% This function stores the gaussian points and the corresponding weights
% for a 1D element.
% 
%   Input:
%       nGauss      : Number of Gauss points used for integration.
%
%   Output:
%       gWeight     : The matrix storing the weights of the Gaussian
%                     integration.
%       gCord       : The vector containing the Gaussian coordinates for
%                     integration.
%   
%
function [gWeight, gCord] = getGaussCordinatesOn1DElement( nGauss )



if (nGauss == 0)
              gWeight = [0];
              gCord(:,:,1) = [0];
elseif(nGauss == 2)
              gWeight = [1, 1];
              gCord(:,:,1) = [-0.57735];
              gCord(:,:,2) = [0.57735];
elseif(nGauss == 3)
              gWeight = [-0.555555, 0.888889, 0.5555555];
              gCord(:,:,1) = [-0.774597];
              gCord(:,:,2) = [0];
              gCord(:,:,3) = [0.774597];

end

