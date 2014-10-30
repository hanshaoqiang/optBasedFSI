%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_referenceConfigurationIncompressibleFlow2D(p,q,U,V,CP,rb,prb,F)
%% Function documentation
%
% Plots the geometry together with the boundary conditions and loads
% arrows for the case of a 2-Dimensional isogeometric incompressible flow.
%
%   Input : 
%     p,q : The polynomial degrees in u-,v-directions
%     U,V : The knot vectors in u-,v- directions
%      CP : The set of the Control points and weights
%      rb : The global numbering of the homogeneous boundary conditions
%     prb : The global numbering of the inhomogeneous boundary conditions
%       F : The force vector
%
%  Output : 
%           graphics
%
% Function layout :
%
%% Function main body

% Create arrays needed for the plotting of the surface, the supports as
% well as the force arrows
gridu = 49;
gridv = 49;
[Xp,Yp,Zp] = createBSplineSurface(p,q,U,V,CP,gridu,gridv);
[xs,ys,zs] = createSupportsForIncompressibleFlow2D(CP,rb);
[pxs,pys,pzs] = createSupportsForIncompressibleFlow2D(CP,prb);
[xf,yf,zf] = createForceArrowsForIncompressibleFlow2D(CP,F);

% Assign an index to the figure
% figure(graph.index)

% geometry
surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none');
hold on;
 
% element edges
isDeformed = 0;
createEdges(p,q,U,V,CP,isDeformed,gridu,gridv)
  
% Plot the homogeneous Dirichlet boundary conditions
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end

% Plot the inhomogeneous Dirichlet boundary conditions
if norm(prb)~=0
    for k =1:length(pxs(:,1))
        plot3(pxs(k,:),pys(k,:),pzs(k,:),'Linewidth',2,'Color','red');
    end
end
    
% load arrows
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end
  
% control points and polygon
createControlPolygon(CP)

axis equal;
% view(0,240);
view(2);
camlight left; lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
% xlim([-4 12]);
% ylim([0 1]);
% zlim([-8 1]);
%grid off;
% set(gca,'xtick',[-100 100]);
% set(gca,'ytick',[-100 100]);
% set(gca,'ztick',[-100 100]);
% axis off;
 
% hold off;

end