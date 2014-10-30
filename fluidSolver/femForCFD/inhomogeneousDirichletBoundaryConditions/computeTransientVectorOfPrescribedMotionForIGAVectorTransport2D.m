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
function vrb = computeTransientVectorOfPrescribedMotionForIGAVectorTransport2D(CP,prb,prescribedMotionX,prescribedMotionY,t)
%% Function documentation
%
%             Input :
%               prb :  The vector containing the global numbering of the DoFs
%                     where inhomogeneous Dirichlet boundary conditions are
%                     applied
% prescribedMotionX : The prescibed motion on the inhomogeneous boundary
%                     conditions for the scalar X
% prescribedMotionY : The prescibed motion on the inhomogeneous boundary
%                     conditions for the scalar Y
%                 t : The current simulation time
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the DoFs where inhomogeneous Dirichlet boundary conditions are applied
%
%% Function main body

%% 0. Read input

% Total number of Control Points in  u-direction
nu = length(CP(:,1,1));

% Initialize output
vrb = zeros(1,length(prb));

%% 1. Loop over all the DoFs where inhomogeneous Dirichlet boundary conditions are applied
for k=1:length(prb)
    % Get the corresponding Control Point number p and indices CP(i,j)
    h = prb(k)/2;
    p = ceil(h);
    j = ceil(p/nu);
    i = p-(j-1)*nu;
    
    % Check in which direction the prescribed DoF is (i.e. to which scalar 
    % is that related)
    
    % If prb is odd is the scalar X
    if p~=h
        if isnumeric(prescribedMotionX)
            % Assign the constant value to the prescribed motion vector
            vrb(k) = prescribedMotionX;
        else
            % Get the Cartesian coordinates of the related Control Point
            x = CP(i,j,1);
            y = CP(i,j,2);
            z = CP(i,j,3);
            
            % Assign the functional value to the prescribed motion vector
            vrb(k) = prescribedMotionX(x,y,z,t);
        end
    % If prb is even is the scalar Y
    else
        if isnumeric(prescribedMotionY)
            % Assign the constant value to the prescribed motion vector
            vrb(k) = prescribedMotionY;
        else
            % Get the Cartesian coordinates of the related Control Point
            x = CP(i,j,1);
            y = CP(i,j,2);
            z = CP(i,j,3);
            
            % Assign the functional value to the prescribed motion vector
            vrb(k) = prescribedMotionY(x,y,z,t);
        end
    end 
end

end

