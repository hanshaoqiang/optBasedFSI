function [Global_Stiffness_Matrix, Global_Mass_Matrix] = ComputeGlobalStiffnessMatrix(diffusionCoeff, mesh)
%   This functions computes the global stiffness matrix for the grid
%   Also calculates the elemental stiffness matrices and assembles them.
%
%   nodes         - Is the matrix containing the Nodal Coordinates
%   elements      - Is the matrix containing the Node numbers of the elements
%

% Computing number of elements and nodes
Number_of_Nodes = max(length(mesh.nodes));
Number_of_Elements = max(length(mesh.elements));

% Declaring the Golobal Stiffness Matrix and local Stiffness matrix
Global_Stiffness_Matrix = zeros(2*Number_of_Nodes, 2*Number_of_Nodes);
Global_Mass_Matrix = zeros(2*Number_of_Nodes, 2*Number_of_Nodes);

% Looping over all the elements to find the Elemental stiffness matrix and
% then assembling it to the Global Stiffness Matrix.
for Elem_Number = 1:Number_of_Elements
    
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
    
    
    %% Computing the stiffness matrix associated with diffusion term
    Element_Diffusion_Matrix = ComputeElementDiffusionMatrix(diffusionCoeff, mesh,Elem_Number);
    Element_Mass_Matrix = ComputeElementMassMatrix(mesh,Elem_Number);
    
    %% Final element matrix
    Element_Matrix = Element_Diffusion_Matrix;
    
   
    %% This process makes it suitable to use Global Stiffness Matrix directly in
    % accumilating the Global Stiffnessmatrix. This is nothing but
    % accumilating the Element stiffness matrix in to Global Stiffness
    % Matrix.
    %
    % Element freedom table
    EFT = zeros(1,6);
    
    % Assign the entried of the element freedom table recursively
    for j=1:3
        EFT(1,2*j-1) = 2*element_vertices(j)-1;
        EFT(1,2*j) = 2*element_vertices(j);
    end
    
    Global_Stiffness_Matrix(EFT,EFT) = Global_Stiffness_Matrix(EFT, EFT) + Element_Matrix;
    Global_Mass_Matrix(EFT,EFT) = Global_Mass_Matrix(EFT, EFT) + Element_Mass_Matrix;
    
    %End of the Loop over all the elements
end

Global_Stiffness_Matrix = sparse(Global_Stiffness_Matrix);
Global_Mass_Matrix = sparse(Global_Mass_Matrix);
% End of the Function
end