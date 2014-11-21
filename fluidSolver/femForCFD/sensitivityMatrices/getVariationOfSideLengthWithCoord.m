function [dhdX] = getVariationOfSideLengthWithCoord( nodes )
%   This method returns the minimum of all the lengths of the given
%   triangle
%   
%   Input:
%       nodes   -- Matrix containing the nodes of the current tringular
%                   element. Its is a 3x3 matrix.
%       GP      -- Gauss point where the matrices are needed.
%   Output:
%       h      -- minimum side length of the triangle

% Predefining the derivative
DhDp = @(p,pp)([p(1)-pp(1) p(2)-pp(2) p(3) - pp(3)]);
DhDpp = @(p,pp)([-p(1)+pp(1) -p(2)+pp(2) -p(3)+pp(3)]);

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


% Side 1 p1 to p2 
l1 = sqrt( (x1-x2)^2 + (y1-y2)^2 );
% Side 2 p2 to p3
l2 = sqrt( (x2-x3)^2 + (y2-y3)^2 );
% Side 3 p3 to p1 
l3 = sqrt( (x1-x3)^2 + (y1-y3)^2 );

[h,n] = min([l1,l2,l3]);

if(n == 1)
    dhdp = DhDp([x1 y1 0], [x2 y2 0]);
    dhdpp = DhDpp([x1 y1 0], [x2 y2 0]);
    dhdX = (1.0/h)*[dhdp(1) dhdp(2) dhdpp(1) dhdpp(2) 0 0]';
elseif(n == 2)
    dhdp = DhDp([x2 y2 0], [x3 y3 0]);
    dhdpp = DhDpp([x2 y2 0], [x3 y3 0]);
    dhdX = (1.0/h)*[0 0 dhdp(1) dhdp(2) dhdpp(1) dhdpp(2)]';
elseif(n == 3)
    dhdp = DhDp([x1 y1 0], [x3 y3 0]);
    dhdpp = DhDpp([x1 y1 0], [x3 y3 0]);
    dhdX = (1.0/h)*[dhdp(1) dhdp(2) 0 0 dhdpp(1) dhdpp(2)]';
end

% End of the Function        
end

