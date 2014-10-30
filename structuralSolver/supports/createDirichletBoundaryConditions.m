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
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rb = createDirichletBoundaryConditions(p,U,q,V,CP,boundaryConditions,mesh)
%% Function documentation
%
% Computes the boundary conditions vector rb given the mesh boundary nodes
%
%              Input :
%                p,q : NURBS polynomial degree
%                U,V : Knot vectors in u,v-direction
%                 CP : set of control point coordinates and weights
% boundaryConditions : contains information on the BC's:
%                             increment : Increment on the search algorithm
%                      tolerance_points : Tolerance on the search algorithm
%                       dirichlet.dofs :
%                                    1 : Prescribe DoF ux
%                                    2 : Prescribe DoF uy   
%                                    3 : Prescribe DoF ux and uy   
%               mesh : Contains nodes, elements boundary keypoints and 
%                      boundary nodes of the mesh
%
%             Output :
%                 rb : Vector 1xn containing the numbering of the
%                      prescribed DoF's
%
%% Function main body

% Number of Control Points in u-direction
nu = length(CP(:,1,1));

% Number of Control Points in v-direction
nv = length(CP(1,:,1));

% Initialize counter for the Dirichlet boundary conditions
rb_counter = 1;

% Loop over all the boundary segments with Dirichlet BC's
for i=1:boundaryConditions.dirichlet.number_of_segments
    
    % Get the parametric domain of the segment i
    segmentu = boundaryConditions.dirichlet.segmentu(i,:);
    segmentv = boundaryConditions.dirichlet.segmentv(i,:);
    
    if isscalar(segmentu)
        % Set up the starting coordinates
        u = segmentu;
        v = segmentv(1);
        
        % Fixed coordinate
        is_on_u = 0;
        
        % Running coordinate
        s = v;
        S = segmentv;
        
    elseif isscalar(segmentv);
        % Set up the starting coordinates
        u = segmentu(1);
        v = segmentv;
        
        % Fixed coordinate
        is_on_u = 1;
        
        % Running coordinate
        s = u;
        S = segmentu;
    end
    
    % Initialize Dirichlet point counter
    dirichlet_point_counter = 1;
    
    % Loop over all the test points on the Dirichlet domain
    while s<=S(length(S))+boundaryConditions.tolerance_search
        
        % Find the span index in u-direction
        spanu = findKnotSpan(u,U,nu);

        % Find the span index in u-direction
        spanv = findKnotSpan(v,V,nv);
        
        % Get the point on Dirichlet domain
        point_on_dirichlet(dirichlet_point_counter,:) = computePointCartesianCoordinatesOnBSplineSurface(p,spanu,u,U,q,spanv,v,V,CP);
        
        % update the parametric location
        s = s + boundaryConditions.increment;
        
        if is_on_u
            u = s;
        else
            v = s;
        end
        
        % Update counter
        dirichlet_point_counter = dirichlet_point_counter + 1;
    end
    
    % Find the Dirichlet boundary condition numbering
    for j=1:length(mesh.nodesOnBoundary) 
        for k=1:length(point_on_dirichlet)
            
            % Test if the node on the boundary belongs to the Dirichlet domain under consideration
            node_on_mesh_boundary = mesh.nodes(mesh.nodesOnBoundary(j),:);
            point_on_boundary = point_on_dirichlet(k,:);
            
            % Get their distance in the L2-norm
            L2normDistance = norm(node_on_mesh_boundary-point_on_boundary);
            
            if L2normDistance<=boundaryConditions.tolerance_points
                % Assign the boundary condition for the corresponding DoF
                if boundaryConditions.dirichlet.dofs(i)==1||boundaryConditions.dirichlet.dofs(i)==3
                    % Blocked in x-direction
                    rb(rb_counter) = 2*mesh.nodesOnBoundary(j)-1;
                    
                    % Update Dirichlet BC's counter
                    rb_counter = rb_counter+1;
                end
                if boundaryConditions.dirichlet.dofs(i)==2||boundaryConditions.dirichlet.dofs(i)==3
                    % Blocked in y-direction
                    rb(rb_counter) = 2*mesh.nodesOnBoundary(j);
                    
                    % Update Dirichlet BC's counter
                    rb_counter = rb_counter+1;
                end
                
            end
            
        end        
    end
    
    % Delete dublicate elements from the Dirichlet array
    rb = unique(rb);
    
    % Adjust accordingly the counter
    rb_counter = length(rb)+1;
end


end

