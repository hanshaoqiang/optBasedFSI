function c = cross_product(a,b)
% Function documentation 
%
% returns cross product c = a x b
%
%   Input : 
%     a,b : two vectors
%
%  Output :
%       c : The cross product of a and b
%
%% Function main body

c(1,1) = a(2)*b(3) - a(3)*b(2);
c(2,1) = a(3)*b(1) - a(1)*b(3);
c(3,1) = a(1)*b(2) - a(2)*b(1);

end