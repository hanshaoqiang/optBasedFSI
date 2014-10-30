function [Element_Mass_Matrix] = ComputeElementMassMatrix(grid,Element_Number)
%   This function calculates the element mass matrix for an Element
%   specified throught Element_Number.
%
%   mesh            - Is the mesh
%   Element_Number  - Is the element number for which the Element stiffness matrix should be Computed.
%
Element_Mass_Matrix = zeros(6);
%% Obtaining the Guassian weights and the Shape Function Matrix at those points
[N, GW, cord] = getShapeFunctionsAndWeights(3);

%% Calculating the terms of the matrix using Gaussian Integration Rule
for i = 1:3
    % Divergence of Shape Function Matrix.
      
    Element_Mass_Matrix = Element_Mass_Matrix + (N(:,:,i)'*N(:,:,i))*GW(i);
    
end
   % End of the Function
end
