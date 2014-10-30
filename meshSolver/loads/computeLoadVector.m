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
function FGlobal = computeLoadVector(p,U,q,V,CP,boundaryConditions,materialProperties,mesh,t)
%% Function documentation
%
% Computes the global load vector for the given mesh and the given Neumann
% boundary conditions
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
%                  t : Time of simulation
% materialProperties : The material properties of the structure
%               mesh : Contains nodes, elements boundary keypoints and 
%                      boundary nodes of the mesh
%
%             Output :
%            FGlobal : Global force vector
%
%% Function layout

% The thickness of the plate
thickness = materialProperties.t;

% Total number of elements in the mesh
no_elements = length(mesh.elements);

% Total number of nodes in the mesh
no_nodes = length(mesh.nodes);

% Element number of nodes
no_nodes_element = 3;

% Total number of DoFs in the mesh
no_dofs = 2*no_nodes;

% Element number of DoFs
no_dofs_element = 2*no_nodes_element;

% Number of Control Points in u-direction
nu = length(CP(:,1,1));

% Number of Control Points in v-direction
nv = length(CP(1,:,1));

% Initialize the structure node array
node = struct([]);

% Initialize all nodes to be of type "unloaded"
for i=1:no_nodes
    % For all individual load cases
    for j=1:boundaryConditions.neumann.number_of_segments
        node(i).neumann.load(j).type = 'unloaded';
    end
end

% Loop over all the boundary segments with Dirichlet BC's
for i=1:boundaryConditions.neumann.number_of_segments
    
    % Get the parametric domain of the segment i
    segmentu = boundaryConditions.neumann.segmentu(i,:);
    segmentv = boundaryConditions.neumann.segmentv(i,:);
    
    if isscalar(segmentu) && ~isscalar(segmentv)
        % Set up the starting coordinates
        u = segmentu;
        v = segmentv(1);
        
        % Fixed coordinate
        is_on_u = 0;
        
        % Running coordinate
        s = v;
        S = segmentv;
        
    elseif isscalar(segmentv) && ~isscalar(segmentu)
        % Set up the starting coordinates
        u = segmentu(1);
        v = segmentv;
        
        % Fixed coordinate
        is_on_u = 1;
        
        % Running coordinate
        s = u;
        S = segmentu;
    elseif isscalar(segmentv)  && isscalar(segmentu)
        % Set up the starting coordinates
        u = segmentu;
        v = segmentv;
        
        % Fixed coordinate
        is_on_u = boundaryConditions.neumann.isOnU(i);
        
        % Running coordinate
        if is_on_u
            s = u;
            S = segmentu;
        else
            s = v;
            S = segmentv;
        end
    end
    
    % Initialize Dirichlet point counter
    neumann_point_counter = 1;
    
    % Loop over all the test points on the Dirichlet domain
    while s<=S(length(S))+boundaryConditions.tolerance_search
        
        % Find the span index in u-direction
        spanu = findKnotSpan(u,U,nu);

        % Find the span index in u-direction
        spanv = findKnotSpan(v,V,nv);
        
        % Get the point on Dirichlet domain
        point_on_neumann(neumann_point_counter,:) = computePointCartesianCoordinatesOnBSplineSurface(p,spanu,u,U,q,spanv,v,V,CP);
        
        % update the parametric location
        s = s + boundaryConditions.increment;
        
        if is_on_u
            u = s;
        else
            v = s;
        end
        
        % Update counter
        neumann_point_counter = neumann_point_counter + 1;
    end
    
    % Initialize counter
    l = 1;
    
    % Find the Dirichlet boundary condition numbering
    for j=1:length(mesh.nodesOnBoundary) 
        for k=1:length(point_on_neumann(:,1))
            
            % Test if the node on the boundary belongs to the Dirichlet domain under consideration
            node_on_mesh_boundary = mesh.nodes(mesh.nodesOnBoundary(j),:);
            point_on_boundary = point_on_neumann(k,:);
            
            % Get their distance in the L2-norm
            L2normDistance = norm(node_on_mesh_boundary-point_on_boundary);
            
            if L2normDistance<=boundaryConditions.tolerance_points && ((~isscalar(segmentu) || ~isscalar(segmentv)) || l==1)
                
                node(mesh.nodesOnBoundary(j)).neumann.load = boundaryConditions.neumann.load(i);
                
                % Update counter
                l = l + 1;
            end
        end        
    end
end

% Initialize global load vector
FGlobal = zeros(no_dofs,1);

% Loop over all the elements
for i=1:no_elements

    % Initialize element load vector
    F_element = zeros(no_dofs_element,1);
    
    % Get the current element on the mesh
    element = mesh.elements(i,:);
    
    % Get the nodes of the mesh in a counter-clockwise fashion
    %
    %        i ------------- k
    %          \           /
    %           \         /
    %            \       /
    %             \     /  
    %              \   /
    %               \ /
    %                j
    %
    % The nodes indices
    index_i = element(1);
    index_j = element(2);
    index_k = element(3);
    
    % The coordinates of the nodes
    node_i = mesh.nodes(index_i,:);
    node_j = mesh.nodes(index_j,:);
    node_k = mesh.nodes(index_k,:);
    
    % Get the load types at the nodes
    load_type_at_node_i =  node(index_i).neumann.load.type;
    load_type_at_node_j =  node(index_j).neumann.load.type;
    load_type_at_node_k =  node(index_k).neumann.load.type;
    
    % Check if there is any BC in the current element and compute the element load vector
    if strcmp(load_type_at_node_i,'surfaceLoad')&&strcmp(load_type_at_node_j,'surfaceLoad')
        % Get the length of the boundary ij
        Lij = norm(node_i-node_j);
        
        % Get the value of the load
        load_amplitude_i = node(index_i).neumann.load.value;
        load_amplitude_j = node(index_j).neumann.load.value;

        % Compute the element load vector
        F_element(1:2,1) = 1/2*Lij*thickness*load_amplitude_i;
        F_element(3:4,1) = 1/2*Lij*thickness*load_amplitude_j;
        
    elseif strcmp(load_type_at_node_i,'surfaceLoad')&&strcmp(load_type_at_node_k,'surfaceLoad')
        % Get the length of the boundary ik 
        Lik = norm(node_i-node_k);
        
        % Get the value of the load
        load_amplitude_i = node(index_i).neumann.load.value;
        load_amplitude_k = node(index_k).neumann.load.value;

        % Compute the element load vector
        F_element(1:2,1) = 1/2*Lik*thickness*load_amplitude_i;
        F_element(5:6,1) = 1/2*Lik*thickness*load_amplitude_k;
        
    elseif strcmp(load_type_at_node_j,'surfaceLoad')&&strcmp(load_type_at_node_k,'surfaceLoad')
        % Get the length of the boundary jk
        Ljk = norm(node_j-node_k);
        
        % Get the value of the load
        load_amplitude_j = node(index_j).neumann.load.value;
        load_amplitude_k = node(index_k).neumann.load.value;

        % Compute the element load vector
        F_element(3:4,1) = 1/2*Ljk*thickness*load_amplitude_j;
        F_element(5:6,1) = 1/2*Ljk*thickness*load_amplitude_k;
    elseif strcmp(load_type_at_node_i,'lineLoad') || strcmp(load_type_at_node_j,'lineLoad') || strcmp(load_type_at_node_k,'lineLoad')
        if strcmp(load_type_at_node_i,'lineLoad')
            % Get the value of the load
            load_amplitude_i = node(index_i).neumann.load.value;

            % Compute the element load vector
            FGlobal(2*index_i-1:2*index_i,1) =FGlobal(2*index_i-1:2*index_i,1) + thickness*load_amplitude_i;
        end
        if strcmp(load_type_at_node_j,'lineLoad')
            % Get the value of the load
            load_amplitude_j = node(index_j).neumann.load.value;

            % Compute the element load vector
            FGlobal(2*index_j-1:2*index_j,1) =FGlobal(2*index_j-1:2*index_j,1) + thickness*load_amplitude_j;
        end
        if strcmp(load_type_at_node_k,'lineLoad')
            % Get the value of the load
            load_amplitude_k = node(index_k).neumann.load.value;

            % Compute the element load vector
            FGlobal(2*index_k-1:2*index_k,1) =FGlobal(2*index_k-1:2*index_k,1) + thickness*load_amplitude_k;
        end
    end
    
    % Create the element freedom table for the DoFs
    freedom_table_dofs = zeros(no_dofs_element,1);
    
    for k=1:no_nodes_element
        % Horizontal DoFs
        freedom_table_dofs(2*k-1,1) = 2*element(k)-1;
        
        % Vertical DoFs
        freedom_table_dofs(2*k,1) = 2*element(k);
    end
    
    % Add the contribution to the global load vector via element freedom tables
    FGlobal(freedom_table_dofs,1) = FGlobal(freedom_table_dofs,1) + F_element;
end




%% Chaning the FGlobal with respect to time
FGlobal = FGlobal * 0.5*cos(pi*t);
% FGlobal = FGlobal * sin(pi*10*t);
% FGlobal = FGlobal * t * 0.01;

end

