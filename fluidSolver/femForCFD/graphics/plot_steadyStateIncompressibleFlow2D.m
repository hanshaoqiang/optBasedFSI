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
function index = plot_steadyStateIncompressibleFlow2D(p,q,U,V,CP,rb,prb,F,up,component,graph)
%% Function documentation
%
% Plots two windows: the first one shows the undeformed geometry together
% with the loads and the supports and the second one shows the distribution
% of the chosen resultant over the computational domain.
%
%       Input :
%         p,q : Polynomial degrees
%         U,V : Knot vectors in u,v-direction
%          CP : Control point coordinates and weights
%          rb : Vector containing the global numbering of the homogeneous
%               Dirichlet boundary conditions
%         prb : Vector containing the global numbering of the inhomogeneous
%               Dirichlet boundary conditions
%           F : Force vector
%          up : The nodal solution vector on the Control Points
%   component : Velocities or pressure to be visualized
%       graph : Structure on the graphics
%
%      Output : graphics
%
% Function layout :
%
% 0. Read Input
%
% 1. Draw the first window: Plot the result as a NURBS surface (only if contour plot or streamlines are not requested)
%
% 2. Draw the second window: Color visualization of the results on the first window
%
%% Function main body

%% 0. Read Input

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsIGAIncompressibleFlow2D(CP,up,component);

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of Control Points in u,v-direction
nu = length(CPd(:,1,1));
nv = length(CPd(1,:,1));

% Compute point coordinates for the original geometry
[Xp,Yp,Zp] = createBSplineSurface(p,q,U,V,CP,49,49);

% Compute point coordinates for the deformed geometry, the supports and the
% force arrows
[Xp_def,Yp_def,Zp_def] = createBSplineSurface(p,q,U,V,CPd,49,49);
[xs_def,ys_def,zs_def] = createSupportsForIncompressibleFlow2D(CPd,rb);
[pxs,pys,pzs] = createSupportsForIncompressibleFlow2D(CPd,prb);
[xf,yf,zf] = createForceArrowsForIncompressibleFlow2D(CPd,F);

% Number the current plot
figure(graph.index)

%% 1. Draw the first window: Plot the result as a NURBS surface (only if contour plot or streamlines are not requested)
if 0% component~=5 && component~=6
    subplot(2,1,1);

    % Plot the reference domain together with the element edges
    surf(Xp,Yp,Zp,'FaceColor','no','EdgeColor','none');
    hold on;
    createEdges(p,q,U,V,CP,0,49,49);

    % Plot the resultant on the z-direction over the reference domain together
    % with its edges and the support conditions
    surf(Xp_def,Yp_def,Zp_def,'FaceColor','blue','EdgeColor','none');
    createEdges(p,q,U,V,CPd,0,49,49);
    createControlPolygon(CPd);
    if norm(rb)~=0
        for k =1:length(xs_def(:,1))
            plot3(xs_def(k,:),ys_def(k,:),zs_def(k,:),'Linewidth',2,'Color','black');
        end
    end
    if norm(prb)~=0
        for k =1:length(pxs(:,1))
            plot3(pxs(k,:),pys(k,:),pzs(k,:),'Linewidth',2,'Color','red');
        end
    end

    for k =1:length(xf(:,1))
        plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
        plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
    end

    camlight left; lighting phong;
    axis equal;
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    if component==1
        title('Velocity component u_x');
    elseif component==2
        title('Velocity component u_y');
    elseif component==3
        title('Pressure distribution');
    elseif component==4 || component==5
        title('Velocity magnitude ||u||_2');
    end
    title ('Isogeometric Incompressible Flow in 2D');
    hold off;
end

%% 2. Draw the second window: Color visualization of the results on the first window

% element displacement vectors
upEl = zeros(mu-p-1,mv-q-1,3*(p+1)*(q+1));
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                upEl(i,j,k)   = up(3*(c*nu+b)-2);
                upEl(i,j,k+1) = up(3*(c*nu+b)-1);
                upEl(i,j,k+2) = up(3*(c*nu+b));
                k = k + 3;
            end
        end
    end
end

tol=10e-10;
% counting index of lines
l=1;  
%incremental step for v
if component~=5
    r = (V(mv)-V(1))/49;
else
    r = (V(mv)-V(1))/29;
end
v=V(1);
while v <= V(mv)+tol
    % Find the span in v-direction
    j = findKnotSpan(v,V,nv);
    
    %incremental step for u
    if component~=5
        s = (U(mu)-U(1))/49;
    else
        s = (U(mu)-U(1))/29;
    end
    u = U(1);
    
    % Update the counter in the u-direction
    k = 1;
    while u <= U(mu)+tol
        % Find the span in u-direction
        i = findKnotSpan(u,U,nu);
        
        % Compute the Cartesian coordinates of the current parametric
        % location on the reference NURBS patch
        P(k,l,1:3) = computePointCartesianCoordinatesOnBSplineSurface(p,i,u,U,q,j,v,V,CP);
        
        % Compute the element nodal transport vector
        upActual(:,1) = upEl(i,j,:);
        
        % Compute the actual transport vector at the parametric location
        if component==1 || component==2 || component==3 || component==5 || component==6 
            upVector(k,l,:) = computeNodalVectorIncompressibleFlow2D(i,p,u,U,j,q,v,V,CP,upActual);
        elseif component==4
            uVectorTemp = computeNodalVectorIncompressibleFlow2D(i,p,u,U,j,q,v,V,CP,upActual);
            upVector(k,l) = sqrt(uVectorTemp(1)^2 +uVectorTemp(2)^2);
        end
        
        % Update the counter and the parametic location in u-direction
        k = k + 1;
        u = u + s;
    end
    
    % Update the counter and the parametic location in u-direction
    l = l + 1;
    v = v + r;
end

if 0 %component~=5 && component~=6
    % plot
    subplot(2,1,2);    
end

if component==1 || component==2 || component==3
    surf(P(:,:,1),P(:,:,2),P(:,:,3),upVector(:,:,component));
elseif component==4
    surf(P(:,:,1),P(:,:,2),P(:,:,3),upVector);
elseif component==5
    scale = 4;
%     quiver(P(:,:,1),P(:,:,2),upVector(:,:,1),upVector(:,:,2),'autoscale','off');
    quiver(P(:,:,1),P(:,:,2),upVector(:,:,1),upVector(:,:,2),scale);
elseif component==6
    % Compute the mesh grid
%     [sx,sy] = meshgrid(P(1,1,1),min(min(P(:,:,2))):.1:max(max(P(:,:,2))));
    [sx,sy] = meshgrid(min(min(P(:,:,1))):.01:max(max(P(:,:,1))),min(min(P(:,:,2))):.01:max(max(P(:,:,2))));
    streamline(P(:,:,1),P(:,:,2),-P(:,:,2),P(:,:,1),sx,sy);
    
    % Visualize the streamline
%     streamline(stream2(P(:,:,1),P(:,:,2),uVector(:,:,1),uVector(:,:,2),sx,sy));
    streamline(stream2(P(:,:,1),P(:,:,2),-P(:,:,2),P(:,:,1),sx,sy));
end

if component~=5 && component~=6
    % Graphics options
    shading interp;
    colormap('default');
    
    % On the color bar
    colorbar;
    hold on;
end

% invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);

% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);

% Create the element edges on the reference NURBS patch
if component~=5 && component~=6
    createEdges(p,q,U,V,CP,0,49,49);
end

view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if component==1
    title('Velocity component u_x');
elseif component==2
    title('Velocity component u_y');
elseif component==3
    title('Pressure distribution');
elseif component==4
    title('Velocity magnitude ||u||_2');
elseif component==5
    title('Contours of the velocity field u');
end

hold off;
index = graph.index + 1;
end