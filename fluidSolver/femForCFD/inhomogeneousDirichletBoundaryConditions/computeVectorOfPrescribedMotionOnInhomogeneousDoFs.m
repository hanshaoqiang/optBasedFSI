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
function [prb,vrb] = computeVectorOfPrescribedMotionOnInhomogeneousDoFs(prb,vrb,u,v,direction,prescribedMotion,CP)
%% Function documentation
%
%               Input :
%                 prb : The outdated global numbering of the DoFs with 
%                       inhomogeneous Dirichlet boundary conditions
%                 vrb : The outdated vector of the prescribed motion on the
%                       DoFs with inhomogeneous Dirichlet boundary
%                       conditions
%                 u,v : region to be fixed (e.g. u = [0 1], v = [1 1])
%           direction : Direction of the prescribed motion
%    prescribedMotion : Prescribed motion in the chosen direction
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
% 2. Loop over the supports and delete double entries if any
%
%% Function main body

%% 0. Read input

% counters
prbCounter = length(prb) + 1;
vrbCounter = length(vrb) + 1;

% Check if the given input is consistent
if (prbCounter~=vrbCounter); error('The given vectors do not match'); end
    
% number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

%% 1. Iterate and add new supports preserving the old ones
for j = v(1)*(nv-1)+1:v(2)*(nv-1)+1
    for i = u(1)*(nu-1)+1:u(2)*(nu-1)+1
        % Get the DoF to be prescribed
        prb(prbCounter) = 2*((j-1)*nu + i-1) + direction;
        
        % Round to nearest integer
        prb(prbCounter) = round(prb(prbCounter));
        
        % Add the prescribed to the DoF value
        if isscalar(prescribedMotion)
            vrb(prbCounter) = prescribedMotion;
        end
        
        % Update counter
        prbCounter = prbCounter + 1;
    end
end

%% sort rb and delete double entries 

% Sort the values
[prb,indexAray] = sort(prb);
vrb = vrb(indexAray);

% Initialize counter
i=1;
% Loop over the supports and delete double entries if any
while i < length(prb)
    if prb(i)==prb(i+1)
        % Delete the double entry in the prb array
        prb(i+1) = [];  
        
        % Delete the double entry in the vrb array but add its value to the
        % prescribed motion
        vrb(i) = vrb(i) + vrb(i+1);
        vrb(i+1) = [];
        
        % Decrease counter
        i=i-1;  
    end
    % Update counter
    i=i+1;
end
    
end

