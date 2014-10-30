
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
function index = plot_colorDistributionForNSEquations2D(mesh,up,graph,component)
%% Function documentation
%
% Plots the distribution of the scalar over the computational domain for
% the 2D convection-diffusion equation using a color for the distribution
% of the convected scalar quantity
%
%   Input :
%    mesh : Structure containing the nodes and the elements of the
%           triangular mesh
%  scalar : The nodal solution vector to the 2D scalar transport equation
%           (2D advection-diffussion)
%   graph : On the graphics
%
%  Output :
%   index : The index of the current graph
%
% Function layout :
%
% 1. Graphical output of the mesh
%
% 2. Loop over all elements and plot their basis functions
%
%% Function main body

%% 1. Graphical output of the mesh
figure(graph.index)
figure(graph.index+1)

% Plot the mesh edges
plot_EdgesTriangularMesh2D(mesh)

% hold on until to plot all basis functions on the mesh
hold on;

%% 2. Loop over all elements and plot their basis functions

% Grid at each triangular element
grid_xi = 2;
grid_eta = 2;

% step size to both directions
step_xi = 1/(grid_xi+1);
step_eta = 1/(grid_eta+1);

for c_element=1:size(mesh.elements(:,1))

    % Compute the basis functions evaluated at the grid points
    P = zeros(grid_xi+2,grid_eta+2,3);
    scalarOverDomain = zeros(grid_xi+2,grid_eta+2);

    % Initialize local coordinates
    xi = 0;
    eta = 0;
    
    % get element from the mesh
    element = mesh.elements(c_element,:);
    
    % get coordinates of the element vertices
    node_i = element(1,1);
    node_j = element(1,2);
    node_k = element(1,3);
    
    % The vertices of the current triangle
    Pi = mesh.nodes(node_i,:);
    Pj = mesh.nodes(node_j,:);
    Pk = mesh.nodes(node_k,:);
    
    
    % Element Freedom Table
    EFT = zeros(1,9);
    % Assign the values to the EFT
    for j=1:3
       EFT(1,3*j-2) = 3*element(1,j)-2;
       EFT(1,3*j-1) = 3*element(1,j)-1;
       EFT(1,3*j) = 3*element(1,j);
    end
    
    if(component==4)
        % Get the element scalar vector that is the magnitude of velocity
        for j=1:3
            u(j) = up( EFT(1,3*j-2) );
            v(j) = up( EFT(1,3*j-1) );
        end
        scalarElement = sqrt(u.*u + v.*v);
    end
    
    if(component==5)
        for j=1:3
            scalarElement(j) = up( EFT(1,3*j) );
        end
    end

    % The moving vertices
    Pi_eta = Pi;
    Pj_eta = Pj;
    
    % Loop over all the sampling points
    for j=1:grid_eta+2
        for i=1:grid_xi+2
            % Get the point in the interior of the line defined from Pi and Pj
            P(i,j,:) = xi*Pi_eta+(1-xi)*Pj_eta;
        
            % Evaluate the CST basis functions at P
            N = computeBasisFunctionsTriangle2D(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
            
            % Evaluate the scalar field value on P by looping over all the
            % products of the basis functions with the nodal solution
            for k = 1:3
                scalarOverDomain(i,j) = scalarOverDomain(i,j) + N(k) * scalarElement(k);
            end
                  
            % Update local coordinate by the step size
            xi = xi + step_xi;
        end
        % Reset xi local coordinate
        xi =0;
    
    	% Update eta local coordinate
        eta = eta + step_eta;
    
        % Update the moving vertices over the triangle edges
        Pi_eta = eta*Pk+(1-eta)*Pi;
        Pj_eta = eta*Pk+(1-eta)*Pj;
    end


    % Graphical output  of the basis functions
    surf(P(:,:,1),P(:,:,2),P(:,:,3),scalarOverDomain(:,:,1),'EdgeColor','none');
    hold on;
end

% Graphics options
colormap('default');
colorbar;
hold on;

% graph properties
view (2);
axis equal;
axis on;
grid on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Concetration of the Transported Scalar Field over the Computation Domain');

hold off;

% Update graph index
graph.index = graph.index + 1;
index = graph.index;
end

