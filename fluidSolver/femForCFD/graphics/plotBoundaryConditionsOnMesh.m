%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               National Technical University of Athens                   %
%               _______________________________________                   %
%                                                                         %
%                 Laboratoty of Testing and Materials                     %
%                 ___________________________________                     %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Efstathios E. Theotokoglou, Professor NTUA                            %
%   (stathis@central.ntua.gr)                                             %
%                                                                         %
%   Andreas Apostolatos, Research Associate TUM                           %
%   (andreas.apostolatos@tum.de)                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plotBoundaryConditionsOnMesh(mesh,rb,F_global,graph)
%% Function documentation
%
% Plots the boundary conditions (Dirichlet/Neumann) onto a given mesh
%
%    Input :
%     mesh : Nodes and elements of the mesh
%       rb : Vector containing the the DoFs (via their global numbering)
%            which are prescribed 
% F_global : Global load vector
%    graph : On the graphics
%
%   Output :
%    index : The index of the current graph
%
%% Function main body

% Create the supports
[xs,ys,zs] = createSupports(mesh.nodes,rb);

%supports
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end

% Hold on the graphics
hold on;

% Create the force arrows
[xf,yf,zf] = createForceArrows(mesh.nodes,F_global);

% load arrows  
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% title of the graph
title('The initial mesh of the reference configuration');
axis on;

% Update the graph index
index = graph.index + 1;

end

