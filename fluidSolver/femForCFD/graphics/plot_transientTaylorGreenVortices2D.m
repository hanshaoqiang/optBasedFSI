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
function plot_transientTaylorGreenVortices2D(p,U,q,V,CP,parameters,transientAnalysis,component,graph)
%% Function documentation
%
% Plots the Taylor-Green vortices over the time for a 2D incompressible
% flow.
%
%             Input :
%               p,q : Polynomial degrees of the 2D NURBS patch
%               U,V : The knot vectors of the NURBS patch
%                CP : The set of Control Point coordinates and weights for the NURBS
%                     patch
%        parameters : The parameters of the flow (density, viscosity)
% transientAnalysis : Transient analysis parameters : 
%                       Tstart : Start time of the simulation
%                    	  Tend : End time of the simulation
%                           nT : Number of time steps
%                        	dt : Time step (automatically computed)
%         component : Velocities or pressure to be visualized
%             graph : Structure on the graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the time steps
%
%    1i. Preamble of the time stepping iterations
%
%        1iv.1. Get the current Cartesian location
%
%        1iv.2. Compute the resultant on the evaluation point
%
%        1iv.3. Update the Cartesian loaction along the x-direction
%
%    1v. Visualize the selected resultant over the domain
%
%   1vi. Update the time of the simulation
%
% 2. Appendix
%
%% Function main body
fprintf('______________________________________________\n');
fprintf('\n');
fprintf('Analytical transient simulation of the Taylor-\n');
fprintf('Green vortices for the Stokes problem has been\n');
fprintf('initiated.\n');
fprintf('______________________________________________\n');
fprintf('\n');

%% 0. Read input

% Number of Control Points
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Number of knots
mu = length(U);
mv = length(V);

% Check input
checkBSplineInput2D(p,mu,nu,q,mv,nv);

% Initialize time
t = transientAnalysis.Tstart;

% Get the vertices of the domain

% Lower left corner
u = 0;
v = 0;
uSpan = findKnotSpan(u,U,nu);
vSpan = findKnotSpan(v,V,nv);
x0 = computePointCartesianCoordinatesOnBSplineSurface(p,uSpan,u,U,q,vSpan,v,V,CP);

% Lower right corner
u = 1;
v = 0;
uSpan = findKnotSpan(u,U,nu);
vSpan = findKnotSpan(v,V,nv);
x1 = computePointCartesianCoordinatesOnBSplineSurface(p,uSpan,u,U,q,vSpan,v,V,CP);

% Upper left corner
u = 0;
v = 1;
uSpan = findKnotSpan(u,U,nu);
vSpan = findKnotSpan(v,V,nv);
x2 = computePointCartesianCoordinatesOnBSplineSurface(p,uSpan,u,U,q,vSpan,v,V,CP);

% Number of evaluation points
nX = 49;
nY = 49;

% Length step at x- and y- direction
dx = (x1-x0)/nX;
dy = (x2-x0)/nY;

% Initialize string
reverseStr = '';

%% 1. Loop over all the time steps
for i=1:transientAnalysis.nT
    %% 1i. Preamble of the time stepping iterations
    
    % Print message on the progress of the load steps
    msg = sprintf('\t Time step %d/%d at real time %d seconds \n',i,transientAnalysis.nT,t);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    %% 1ii. Initialize the arrays
    
    % Initialize the field arrays
    if component~=5
        resultant = zeros(nX,nY);
    else
        resultant = zeros(nX,nY,2);
    end
    
    % Initialize the Cartesian components
    PCartesian = zeros(nX,nY,3);
    
    %% 1iii. Initialize the Coordinates on the Cartesian space
    x = x0;
    
    %% 1iv. Loop over all the evaluation points
    for j=1:nX+1
        for k=1:nY+1
            %% 1iv.1. Get the current Cartesian location
            PCartesian(j,k,1:3) = x;
            
            %% 1iv.2. Compute the resultant on the evaluation point
            if component==1
                resultant(j,k) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
            elseif component==2
                resultant(j,k) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
            elseif component==3
                resultant(j,k) = -.25*(cos(2*x(1)) + cos(2*x(2)))*exp(-4*t*parameters.nue);
            elseif component==4
                uX = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
                uY = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
                resultant(j,k) = sqrt(uX^2 + uY^2);
            elseif component==5
                resultant(j,k,1) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
                resultant(j,k,2) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
            end
            
            %% 1iv.3. Update the Cartesian loaction along the x-direction
            x = x + dx;
        end
        %% 1iv.3.. Update the Cartesian loaction along the y-direction
        x = x0;
        x = x + j*dy;
    end
    
    %% 1v. Visualize the selected resultant over the domain
    
    % Initialize figure handle
    figure(graph.index);

    % Plot the resultant over the domain
    if component==1 || component==2 || component==3 || component==4
        surf(PCartesian(:,:,1),PCartesian(:,:,2),PCartesian(:,:,3),resultant(:,:),'EdgeColor','none');
    elseif component==5
        quiver(PCartesian(:,:,1),PCartesian(:,:,2),resultant(:,:,1),resultant(:,:,2),'autoscale','off');
    end
    
    if component~=5 && component~=6
        % Graphics options
        shading interp;
        colormap('default');

        % On the color bar
        colorbar;
        hold on;
    end
    
    % Hold on the graph
    hold on;
    
    % Plot the edges of the elements
    if component~=5 && component~=6
        createEdges(p,q,U,V,CP,0,49,49);
    end

    % On the graph properties
    view(2);
    axis equal;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    if component==1
        title('Velocity component u_x');
    elseif component==2
        title('Velocity component u_y');
    elseif component==3
        title('Pressure distribution');
    elseif component==4
        title('Velocity magnitude ||u||_2');
    elseif component==5
        title('Contours of the velocity field u');
    end
    
    % Hold off the graph
    hold off;
    
    %% 1vi. Update the time of the simulation
    t = t + transientAnalysis.dt;
    
end

% Close the graph
close(graph.index);

%% 2. Appendix
fprintf('\n');
fprintf('__________Transient Simulation Ended__________\n');

end

