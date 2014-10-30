function [ x, y ] = getGaussPointCoordinates( GP, nodes)
%   This method contains routene for calculating the derivatives of Nodal
%   basis on a particular element.
%   
%   Input:
%       nodes   -- Matrix containing the coordinates of nodes of the current tringular
%                   element. Its is a 3x3 matrix.
%       GP      -- Gauss point where the matrices are needed.
%   Output:
%       N       -- Shape function Matrix. (3x9)
%       dN      -- First derivatives of the Shape Functions.(9x9)
%       dDotN   -- Matrix Divergence of N. (1x9)
%       laplaceN-- Laplace of the Shape Functions.[This will be
%                   zero for linear Triangles] (3x9)

% Obtaining the x and y coordinates of all the nodes.
% NODE 1
x1 = nodes(1,1);
y1 = nodes(1,2);
% NODE 2
x2 = nodes(2,1);
y2 = nodes(2,2);
% NODE 3
x3 = nodes(3,1);
y3 = nodes(3,2);

% Computing area
area = 0.5*abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );

%% Calculating all the derivatives of the Shape functions
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


%% Finding the values of the shape functions at gauss points
% Since the Gauss points are already in the parametric coordinates we can
% directly use them. This is a 3x9 matrix
N = [GP(1) 0 0 GP(2) 0 0 GP(3) 0 0;
     0 GP(1) 0 0 GP(2) 0 0 GP(3) 0;
     0 0 0     0 0 0     0 0 0];


x = GP(1)*x1 + GP(2)*x2 + GP(3)*x3;
y = GP(1)*y1 + GP(2)*y2 + GP(3)*y3;

        

% End of the Function        
end

