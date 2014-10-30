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
function index = plot_transientVectorTransportProblem2D(p,q,U,V,CP,rb,F,transportVector,component,graph)
%% Function documentation
%
% Plots the distribution of the chosen scalar quantity over the
% computational domain for a transient analysis of the 2D vector
% convection-diffusion-reaction equation.
%
%       Input :
%         p,q : Polynomial degrees
%         U,V : Knot vectors in u,v-direction
%          CP : Control point coordinates and weights
%          rb : Vector containing the global numbering of the homogeneous
%               Dirichlet boundary conditions
%           F : Force vector
%           u : The scalar field on the Control Points
%   component : Scalar component to be visualized
%       graph : Structure on the graphics
%
%      Output : graphics
%
% Function layout :
%
% 0. Read Input
%
% 1. Draw the first window: Plot the result as a NURBS surface
%
% 2. Draw the second window: Color visualization of the results on the first window
%
%% Function main body

%% 0. Read Input

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsIGAVectorTransport2D(CP,transportVector,component);

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
[xs_def,ys_def,zs_def] = createSupports2D(CPd,rb);
[xf,yf,zf] = createForceArrows2D(CP,F);

% Number the current plot
figure(graph.index)

%% Plot the result as a NURBS surface
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
for k =1:length(xs_def(:,1))
    plot3(xs_def(k,:),ys_def(k,:),zs_def(k,:),'Linewidth',2,'Color','black');
end
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

camlight left; lighting phong;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if component == 1
    zlabel('u_x','FontSize',14);
elseif component == 2
    zlabel('u_y','FontSize',14);
elseif component == 3
    zlabel('||u||_2','FontSize',14);
end
title ('Isogeometric Vector Transport in 2D');
hold off;

%% 2. Draw the second window: Color visualization of the results on the first window

% element displacement vectors
uEl = zeros(mu-p-1,mv-q-1,2*(p+1)*(q+1));
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                uEl(i,j,k)   = transportVector(2*(c*nu+b)-1);
                uEl(i,j,k+1) = transportVector(2*(c*nu+b));
                k=k+2;
            end
        end
    end
end

tol=10e-10;
% counting index of lines
l=1;  
%incremental step for v
r=(V(mv)-V(1))/49;    
v=V(1);
while v <= V(mv)+tol
    % Find the span in v-direction
    j = findKnotSpan(v,V,nv);
    
    %incremental step for u
    s = (U(mu)-U(1))/49;
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
        uActual(:,1) = uEl(i,j,:);
        
        % Compute the actual transport vector at the parametric location
        if component==1 || component==2
            uVector(k,l,:) = computeVectorTransport2D(i,p,u,U,j,q,v,V,CP,uActual);
        elseif component==3
            uVectorTemp = computeVectorTransport2D(i,p,u,U,j,q,v,V,CP,uActual);
            uVector(k,l) = sqrt(uVectorTemp(1)^2 +uVectorTemp(2)^2);
        end 
        
        % Update the counter and the parametic location in u-direction
        k=k+1;
        u=u+s;
    end
    
    % Update the counter and the parametic location in u-direction
    l=l+1;
    v=v+r;
end

% plot
subplot(2,1,2);

if component==1 || component==2
    surf(P(:,:,1),P(:,:,2),P(:,:,3),uVector(:,:,component));
elseif component==3
    surf(P(:,:,1),P(:,:,2),P(:,:,3),uVector);   
end

% Graphics options
shading interp;
colormap('default');

% invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);

% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);

colorbar;
hold on;

% Create the element edges on the reference NURBS patch
createEdges(p,q,U,V,CP,0,49,49);

view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if component == 1
    zlabel('u_x','FontSize',14);
elseif component == 2
    zlabel('u_y','FontSize',14);
elseif component == 3
    zlabel('||u||_2','FontSize',14);
end
if component==1
    title('Vector transport component u_x');
elseif component==2
    title('Vector transport component u_y');
elseif component==3
    title('Vector transport magnitude ||u||_2');
end

hold off;

index = graph.index + 1;
end