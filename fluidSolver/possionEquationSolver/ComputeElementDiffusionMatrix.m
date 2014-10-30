function [Element_Diffusion_Matrix] = ComputeElementDiffusionMatrix(diffusionCoeff, mesh,Elem_Number)
%   This function calculates the element stiffness matrix for an Element
%   specified throught Element_Number.
%
%   nodes         - Is the matrix containing the Nodal Coordinates
%   elements      - Is the matrix containing the Node numbers of the elements
%   Element_Number- Is the element number for which the Element stiffness matrix should be Computed.
%
Element_Diffusion_Matrix = zeros(6);
% Finding the Nodes corresponding to the present element.
element_vertices = mesh.elements(Elem_Number,:);
% Finding the Nodes corresponding to the present element.
element_nodes = [   mesh.nodes(element_vertices(1),:);
    mesh.nodes(element_vertices(2),:);
    mesh.nodes(element_vertices(3),:) ];

% NODE 1
x1 = element_nodes(1,1);
y1 = element_nodes(1,2);
% NODE 2
x2 = element_nodes(2,1);
y2 = element_nodes(2,2);
% NODE 3
x3 = element_nodes(3,1);
y3 = element_nodes(3,2);

area = 0.5*abs(x1*(y2-y3)+ x2*(y3-y1)+x3*(y1-y2));

%% #############
% All the formulation below is based on the equations/literature developed
% by me based on the references especially on SC1 Exam 2012. For Basis
% Functions and Triangle mapping PDF in the Reference Folder.
%
% NOTE : Basis Function used are of Order 1
% #############

% Differential of Basis Function 1
dN1_dx = (y2-y3)/(2*area);
dN1_dy = (x3-x2)/(2*area);

% Differential of Basis Function 2
dN2_dx = (y3-y1)/(2*area);
dN2_dy = (x1-x3)/(2*area);

% Differential of Basis Function 3
dN3_dx = (y1-y2)/(2*area);
dN3_dy = (x2-x1)/(2*area);

%% Obtaining the Guassian weights and the Shape Function Matrix at those points
[N, GW, cord] = getShapeFunctionsAndWeights(3);

%% Calculating the terms of the matrix using Gaussian Integration Rule
for i = 1:3
    % Divergence of Shape Function Matrix.
    B=[dN1_dx,0,dN2_dx,0,dN3_dx,0;
        0,dN1_dx,0,dN2_dx,0,dN3_dx;
        dN1_dy,0,dN2_dy,0,dN3_dy,0;
        0,dN1_dy,0,dN2_dy,0,dN3_dy];
    
    Element_Diffusion_Matrix = Element_Diffusion_Matrix + B'*B.*GW(i);
    
end
Element_Diffusion_Matrix = diffusionCoeff.* Element_Diffusion_Matrix;

% End of the Function
end
