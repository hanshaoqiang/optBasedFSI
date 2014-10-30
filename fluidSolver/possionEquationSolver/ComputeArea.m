function[area] = ComputeArea(grid,Elem_Number)
% This function computes the area of a given Element
%   nodes       - Is matrix containing the nodal Coordinates.
%   elements    - Is matrix containing the nodes of all the elements.
%   Elem_Number - Is the number of element for which area is to be
%                 calculated.

% Obtaining the X and Y coordinates of the nodes for the given elements
node1 = grid.elements(Elem_Number,1); 
node2 = grid.elements(Elem_Number,2);
node3 = grid.elements(Elem_Number,3);
% NODE 1
node1_x = grid.nodes(node1,1);
node1_y = grid.nodes(node1,2);
% NODE 2
node2_x = grid.nodes(node2,1);
node2_y = grid.nodes(node2,2);
% NODE 3
node3_x = grid.nodes(node3,1);
node3_y = grid.nodes(node3,2);


area = 0.5*(node1_x*(node2_y-node3_y)+node2_x*(node3_y-node1_y)+node3_x*(node1_y-node2_y));


if (area < 0)
    area = area *-1;
    display('Area is negative for element ' + Elem_Number);
end

% End of function
end
