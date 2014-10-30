function [Load_Matrix] = ComputeLoadMatrix(physics, Q,grid, tao)
% This function computes and assembles the Right hand side for the system
% of equations that is the LOAD matrix.
%
% nodes         - Is the matrix containing the Nodal Coordinates
% elements      - Is the matrix containing the Node numbers of the elements
% Q(x,y)        - Is the function for the heat generation.
%
% Returns the Assembled load matirx of size (number_of_nodes X 1)
% Extracting the x and y components of the convection velocity
ax = physics.convCoeff.x;
ay = physics.convCoeff.y;

% Computing number of elements and nodes
Number_of_Nodes = length(grid.nodes);
Number_of_Elements = max(size(grid.elements));

% Declaring a zero load matrix
Load_Matrix = zeros(2*Number_of_Nodes,1);

% Looping over All the Elements
for Elem_Number = 1:Number_of_Elements
    
    % Finding the area of the present element
    area = ComputeArea(grid,Elem_Number);
    element_vertices = grid.elements(Elem_Number,:);
    % Finding the Nodes corresponding to the present element.
    node1 = grid.elements(Elem_Number,1);
    node2 = grid.elements(Elem_Number,2);
    node3 = grid.elements(Elem_Number,3);
    
    % Obtaining the x and y coordinates of all the nodes.
    % NODE 1
    x1 = grid.nodes(node1,1);
    y1 = grid.nodes(node1,2);
    % NODE 2
    x2 = grid.nodes(node2,1);
    y2 = grid.nodes(node2,2);
    % NODE 3
    x3 = grid.nodes(node3,1);
    y3 = grid.nodes(node3,2);
    x = [x1, x2, x3];
    y = [y1, y2, y3];
    
    % Differential of Basis Function 1
    dN1_dx = (y2-y3)/(2*area);
    dN1_dy = (x3-x2)/(2*area);
    
    % Differential of Basis Function 2
    dN2_dx = (y3-y1)/(2*area);
    dN2_dy = (x1-x3)/(2*area);
    
    % Differential of Basis Function 3
    dN3_dx = (y1-y2)/(2*area);
    dN3_dy = (x2-x1)/(2*area);
    
    % Element freedom table
    EFT = zeros(1,6);
    %% Obtaining the Guassian weights and the Shape Function Matrix at those points
    [N, GW, cord] = getShapeFunctionsAndWeights(physics.numGP);
    % Jacobian Matrix
    jacobian = ComputeElementJacobianMatrix(grid, Elem_Number);
    % Assign the entried of the element freedom table recursively
    for j=1:3
        EFT(1,2*j-1) = 2*element_vertices(j)-1;
        EFT(1,2*j) = 2*element_vertices(j);
    end
    
    % Performing Gaussian Integration
    for i=1:physics.numGP
        % Coordinates at the present guassian points.
        xGP = sum(x.*cord(:,:,i));
        yGP = sum(y.*cord(:,:,i));
        
        f = [Q.x(xGP,yGP);Q.y(xGP,yGP)];
        n = N(:,:,i)*jacobian;
        k =  n'*f.*GW(i);
        Load_Matrix(EFT) = Load_Matrix(EFT) + k;
        clear f;
        clear k;
    end
    
    %% Calculating the stabilizing term [(a.nabla)N]'*S
    k1 = zeros(6,1);
    for i=1:physics.numGP
        % Coordinates at the present guassian points.
        xGP = sum(x.*cord(:,:,i));
        yGP = sum(y.*cord(:,:,i));
        
        f = [Q.x(xGP,yGP);Q.y(xGP,yGP)];
        B = [(ax*dN1_dx)+(ay*dN1_dy), 0, (ax*dN2_dx)+(ay*dN2_dy), 0, (ax*dN3_dx)+(ay*dN3_dy), 0;
            0, (ax*dN1_dx)+(ay*dN1_dy), 0, (ax*dN2_dx)+(ay*dN2_dy), 0, (ax*dN3_dx)+(ay*dN3_dy)];
        
        k1 = k1 + B'*f.*GW(i);
        clear f;
    end
    Load_Matrix(EFT) = Load_Matrix(EFT) - tao.*k1;
    
    %% Calculating the stabilizing term [N]'*S
    k2 = zeros(6,1);
    for i=1:physics.numGP
        % Coordinates at the present guassian points.
        xGP = sum(x.*cord(:,:,i));
        yGP = sum(y.*cord(:,:,i));
        
        f = [Q.x(xGP,yGP);Q.y(xGP,yGP)];
        k2 =  N(:,:,i)'*f.*GW(i);
        k2 = physics.sigma.*k2;
        
        clear f;
    end
    Load_Matrix(EFT) = Load_Matrix(EFT) + tao.*k2;
    % End of loop over all the elements
end

% End of Funciton
end

