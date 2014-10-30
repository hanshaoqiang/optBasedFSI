%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger   (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                 %
%   Dipl.-Math. Andreas Apostolatos     (andreas.apostolatos@tum.de)       %
%   Aditya Ghantasala                   (aditya.ghantasala@tum.de)        %
%   _______________________________________________________________       %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displacement = getDisplacementAtPoint(mesh, disp, x, y)
%% Documentation of function
%  This method returns the field for Taylor Green Vortices depending on the
%  time of simulation.
%
%   Input :
%       mesh    : The mesh
%       disp      : calculated disp values
%       x       : X Cordinate      
%       y       : Y Cordinate
%   Output :
%       displacement     : displacement field field for taylor green vortices at
%       			requested point x and y
%% Function main body
displacement = zeros(3,1);

elem_num = 0;
%% Loop over all the elements to find out in which element does the point belong to 

for i = 1:length(mesh.elements)
    % Using BAycentric Technique as described in 
    % http://www.blackpawn.com/texts/pointinpoly/
    
    % Obtaining the element nodal coordinates
    element_vertices = mesh.elements(i,:);
    element_nodes = [   mesh.nodes(element_vertices(1),:);
                        mesh.nodes(element_vertices(2),:);
                        mesh.nodes(element_vertices(3),:) ];
    % Computing Vectors
    v0 = element_nodes(3,:) - element_nodes(1,:);
    v1 = element_nodes(2,:) - element_nodes(1,:);
    v2 = [x,y,0] - element_nodes(1,:);
    
    % Computing Dot Products
    dot00 = v0*v0';
    dot01 = v0*v1';
    dot02 = v0*v2';
    dot11 = v1*v1';
    dot12 = v1*v2';
    
    % Computing barycentric Coordinates
    invDenom = 1/(dot00*dot11 - dot01*dot01);
    u = (dot11*dot02 - dot01*dot12)*invDenom;
    v = (dot00*dot12 - dot01*dot02)*invDenom;
    
    
    % Checking if the point is inside the triangle and exiting the loop if
    % it is
    if (u>=0 && v>=0 && (u+v<1))
        elem_num = i;
        break;
    end
    
end

% Calculating the field values at the point using Shape function approach
element_vertices = mesh.elements(elem_num,:);
element_nodes = [   mesh.nodes(element_vertices(1),:);
                        mesh.nodes(element_vertices(2),:);
                        mesh.nodes(element_vertices(3),:) ];


mat = [element_nodes(1,1), element_nodes(2,1), element_nodes(3,1);
        element_nodes(1,2), element_nodes(2,2), element_nodes(3,2);
        1,                  1,                  1 ];
    
    
N = mat \ [x;y;1];

u1 = [disp(2*element_vertices(1)-1); disp(2*element_vertices(1))];
u2 = [disp(2*element_vertices(2)-1); disp(2*element_vertices(2))];
u3 = [disp(2*element_vertices(3)-1); disp(2*element_vertices(3))];
    
displacement = N(1)*u1 + N(2)*u2 + N(3)*u3;
    

% end of the function
end
