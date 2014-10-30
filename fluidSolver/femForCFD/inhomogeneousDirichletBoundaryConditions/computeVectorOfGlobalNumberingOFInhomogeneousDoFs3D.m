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
function prb = computeVectorOfGlobalNumberingOFInhomogeneousDoFs3D(prb,u,v,direction,CP)
%% Function documentation
%
% Returns the vector containing the global numbering of the inhomogeneous
% Dirichlet boundary conditions.
%
%               Input :
%                 prb : The outdated global numbering of the DoFs with 
%                       inhomogeneous Dirichlet boundary conditions
%                 u,v : region to be fixed (e.g. u = [0 1], v = [1 1])
%           direction : Direction of the prescribed motion
%
%              Output : 
%                 vrb : The updated vector of the prescribed motion on the
%                       DoFs with inhomogeneous Dirichlet boundary
%                       conditions
%                 prb : The updated global numbering of the DoFs with
%                       inhomogeneous Dirichlet boundary conditions
%
% Function layout :
%
% 0. Read input
%
% 1. Iterate and add new supports preserving the old ones
%
% 2. Sort rb and delete double entries 
%
%% Function main body

%% 0. Read input

% counters
prbCounter = length(prb) + 1;
    
% number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

%% 1. Iterate and add new supports preserving the old ones
for j = v(1)*(nv-1)+1:v(2)*(nv-1)+1
    for i = u(1)*(nu-1)+1:u(2)*(nu-1)+1
        % Get the DoF to be prescribed
        prb(prbCounter) = 3*((j-1)*nu + i-1) + direction;
        
        % Round to nearest integer
        prb(prbCounter) = round(prb(prbCounter));
        
        % Update counter
        prbCounter = prbCounter + 1;
    end
end

%% 2. Sort rb and delete double entries 
[prb,~] = sort(prb);
prb = unique(prb);
    
end

