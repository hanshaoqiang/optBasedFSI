function [Jacobian] = ComputeElementJacobianMatrix(grid,Element_Number)
%   This function calculates the element Jacobian Matrix for an Element
%   specified throught Element_Number.
%   
%   Element_Number- Is the element number for which the Element stiffness matrix should be Computed.   
%   nodes- is the matrix containing the nodal coordinates
%   elements- is the element connectivity matirx

% Finding the Nodes corresponding to the present element.
node1 = grid.elements(Element_Number,1);
node2 = grid.elements(Element_Number,2);
node3 = grid.elements(Element_Number,3);

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


Jacobian = [1, 0, 1, 0, 1, 0;
            0, 1, 0, 1, 0, 1;
            x1, 0, x2, 0, x3, 0;
            0, x1, 0, x2, 0, x3;
            y1, 0, y2, 0, y3, 0;
            0, y1, 0, y2, 0, y3];


end