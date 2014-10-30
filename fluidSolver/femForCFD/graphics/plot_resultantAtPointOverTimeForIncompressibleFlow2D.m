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
function index = plot_resultantAtPointOverTimeForIncompressibleFlow2D(u,p,U,v,q,V,CP,upHistory,transientAnalysis,component,graph)
%% Function documenation
%
% Returns a 2D plot with the evolution of a resultant throught the
% transient analysis at a specific point corresponding to a 2D
% incompressible flow problem using isogeometric analysis.
%
%                     Input :
%                       u,v : The NURBS parametric coordinates of the
%                             points on which to plot the resultants over 
%                             the time
%                       p,q : The polynomial degrees of the NURBS patch
%                       U,V : The knot vectors of the NURBS patch
%                        CP : The set of Control Point coordinates and
%                             weights for the NURBS patch
%                 upHistory : The history data of the transient analysis by
%                             solving the problem numerically
%         transientAnalysis : Transient analysis parameters : 
%                              Tstart : Start time of the simulation
%                                Tend : End time of the simulation
%                                  nT : Number of time steps
%                                  dt : Time step (automatically computed)
%                component : Velocities or pressure to be visualized
%                    graph : Structure on the graphics
%
%                   Output :
%                            graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Find the knot spans where the point lives in
%
% 2. Compute the Cartesian coordinates of the Point
%
% 3. Loop over all the time steps
%
%    3i. Fill up the time step array accordingly
%
%   3ii. Get the current discrete solution vector
%
%  3iii. Get the actual discrete solution vector
%
%   3iv. Get the actual resultants at the point and at the current time step
%
%    3v. Get the desirable resultant at the point and at the current time step
%
%   3vi. Update the simulation time
%
% 4. Plot the resultant at the given point and over time
%
% 5. Update the figure handle index
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of Control Points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Check input
checkBSplineInput2D(p,mu,nu,q,mv,nv);

% Initialize figure handle
figure(graph.index)

% Initialize the output arrays
resultantAtPointNumerical = zeros(transientAnalysis.nT,1);
timeSteps = zeros(transientAnalysis.nT,1);

% Initialize time
t = transientAnalysis.Tstart;

%% 1. Find the knot spans where the point lives in
uSpan = findKnotSpan(u,U,nu);
vSpan = findKnotSpan(v,V,nv);

%% 2. Compute the Cartesian coordinates of the Point
xCartesian = computePointCartesianCoordinatesOnBSplineSurface(p,uSpan,u,U,q,vSpan,v,V,CP);
x = xCartesian(1);
y = xCartesian(2);

%% 3. Loop over all the time steps
for timeCounter=1:transientAnalysis.nT
    %% 3i. Fill up the time step array accordingly
    timeSteps(timeCounter,1) = t;
    
    %% 3ii. Get the current discrete solution vector
    
    % Distribute the solution vector into the elememts for the current time
    % step :
    
    % element discrete solution vectors
    upEl = zeros(mu-p-1,mv-q-1,3*(p+1)*(q+1));
    for j = (q+1):(mv-q-1)
        for i = (p+1):(mu-p-1)
            % Initialize the counter
            k = 1;
            for c = j-q-1:j-1 
                for b = i-p:i
                    upEl(i,j,k)   = upHistory(3*(c*nu+b)-2,timeCounter);
                    upEl(i,j,k+1) = upHistory(3*(c*nu+b)-1,timeCounter);
                    upEl(i,j,k+2) = upHistory(3*(c*nu+b),timeCounter);

                    % Update the counter
                    k = k + 3;
                end
            end
        end
    end
    
    %% 3iii. Get the actual discrete solution vector
    upActual = upEl(uSpan,vSpan,:);
    upActualVector = zeros(3*(p+1)*(q+1),1);
    for i=1:3*(p+1)*(q+1)
        upActualVector(i) = upActual(1,1,i);
    end
    
    %% 3iv. Get the actual resultants at the point and at the current time step
    upVector = computeNodalVectorIncompressibleFlow2D(uSpan,p,u,U,vSpan,q,v,V,CP,upActualVector);
    
    %% 3v. Get the desirable resultant at the point and at the current time step
    if component==1
        resultantAtPointNumerical(timeCounter,1) = upVector(1);
    elseif component==2
        resultantAtPointNumerical(timeCounter,1) = upVector(2);
    elseif component==3
        resultantAtPointNumerical(timeCounter,1) = upVector(3);
    elseif component==4
        uXNumerical = upVector(1);
        uYNumerical = upVector(2);
        resultantAtPointNumerical(timeCounter,1) = sqrt(uXNumerical^2 + uYNumerical^2);
    end
   
    %% 3vi. Update the simulation time
    t = t + transientAnalysis.dt;
end

%% 4. Plot the resultant at the given point and over time
plot(timeSteps,resultantAtPointNumerical,'b');
legend('Finite Element Approximation');
xlabel('time (seconds)');
if component==1
    yLabelString = 'x-velocity component u_x (m/s)';
elseif component==2
    yLabelString = 'y-velocity component u_y (m/s)';
elseif component==3
    yLabelString = 'pressure p (Pa)';
elseif component==4
    yLabelString = 'velocity magnitude ||u||_2 (m/s)';
end
ylabel(yLabelString);
title(sprintf('Evaluation point X = (%d,%d,%d)',xCartesian(1),xCartesian(2),xCartesian(3)));
    
%% 5. Update the figure handle index
index = graph.index + 1;


end

