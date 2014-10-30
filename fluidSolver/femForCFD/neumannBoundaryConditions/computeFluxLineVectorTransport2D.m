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
function Fl = computeFluxLineVectorTransport2D(Fold,ub,vb,p,q,U,V,CP,int,hFlux,direction)
%% Function documentation
%
% Returns the consistent force vector for the application of a flux over a
% portion of the computational domain's boundary (Neumann boundary) defined
% by parameters ub, vb and which is of magnitude hFlux
% 
%    Input :
%     Fold : The existing force vector
%    ub,vb : Flux application extension (e.g. ub=[0 1], vb=1)
%      p,q : The polynomial degrees of the NURBS description
%      U,V : The knot vectors of the NURBS description
%      int : Structure responsible for the integration
%            automatic : Automatic choice of the Gauss Points fulfiling
%                        stability
%               manual : Manual choice of the Gauss Points for improving
%                        the accuracy
%    fload : constant line load or handle to load function [N/m]
%      dir : direction of h:
%            1=x, 2=y    for h(s)   -> integration over the boundary space
%            3=x, 4=y    for h(x,y) -> integration over the Cartesian space
%            5=parallel, 6=+90  to the edge
%            uv: dir to integrate: 1=u, 2=v. determined from ub,vb
%            xys: 1=integrate ds   2=integrate dx,dy.  determined from dir
%
%   Output :
%       Fl : The updated force vector
%
%% Function main body

% Length of the knot vectors
mu = length(U);
mv = length(V);

% Number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% if the fload parameter is numeric then assign directly its value
if isnumeric(hFlux)==1  
    h = hFlux;  
end

% Issue the Gauss points for the selected integration scheme:
if int.type == 0
   % Default scheme is the full gaussian quadrature element-wise (FGI)
   ngauss(1) = ceil((p+1)/2) ;
   ngauss(2) = ceil((q+1)/2) ;
elseif int.type == 1
    % Manual choice of the gauss points
    ngauss(1) = int.uNgaussLoad;
    ngauss(2) = int.vNgaussLoad;
end

% Find the span in u-direction where to apply the load
i1 = findKnotSpan(ub(1),U,nu);
if (isscalar(ub))
    % If its a scalar one Gauss point sufficies
    u = ub(1);    i2 = i1;
    ugauss = 1;   gwu = 1;
    mapu = 1;     uv = 2;
else
    % If its not a scalar adjust the quadrature respectively
    i2 = findKnotSpan(ub(2),U,nu);
    if ub(2)~=U(mu)    
        i2=i2-1;    
    end
    ugauss = ngauss(1);
end

% Find the span in v-direction where to apply the load
j1 = findKnotSpan(vb(1),V,nv);
if (isscalar(vb))
    % If its a scalar one Gauss point sufficies
    v = vb(1);    j2 = j1;
    vgauss = 1;   gwv = 1;
    mapv = 1;     uv = 1;
else
    % If its not a scalar adjust the quadrature respectively
    j2 = findKnotSpan(vb(2),V,nv);
    if vb(2)~=V(mv)    
        j2=j2-1;    
    end
    vgauss = ngauss(2);
end

% Decide in which coordinate system to integrate
if direction==1||direction==2  
    % Integrate over the curvilinear parametric line
    xys = 1;
elseif direction==3||direction==4||direction==5||direction==6  
    % Integrate over the cartesian coordinate system xy
    xys = 2;    
end

% Integrate the weak form of the load
% Initialize global load vector
F = zeros(nu,nv,2);   

% Initialize element load vector
Fel = zeros(p+1,q+1,2);

% Loop over all elements on which the load is applied
for j = j1:j2
    for i = i1:i2
        % check if we are in a non-zero knot span
        if (U(i+1)~=U(i) && V(j+1)~=V(j))
            % Issue the Gauss points for the correct integration
            [GPu,GWu] = gauss(ngauss(1));
            [GPv,GWv] = gauss(ngauss(2));
        
            % Loop over all Gauss points
            for kv = 1:vgauss
                for ku = 1:ugauss
                    % load extension in u-direction
                    if (isscalar(ub)==0)
                        % map the quadrature point into the knot span in
                        % u-direction
                        u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
                        % compute the respective Jacobian determinant
                        mapu = (U(i+1)-U(i))/2;
                        % issue quadrature weight in u-direction
                        gwu = GWu(ku);
                    end
                    % load extension in v-direction
                    if (isscalar(vb)==0)
                        % map the quadrature point into the knot span in
                        % v-direction
                        v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
                        % compute the respective Jacobian determinant
                        mapv = (V(j+1)-V(j))/2;
                        % issue quadrature weight in v-direction
                        gwv = GWv(kv);
                    end
                    % Compute the product of the Jacobian determinants in
                    % u,v-directions
                    map = mapu*mapv;
                    
                    % Compute the product of the quadrature weights 
                    gw = gwu*gwv;
                    
                    if isnumeric(hFlux) == 0  
                        h = hFlux(p,i,u,U,q,j,v,V,CP);  
                    end
                    
                    % Length of the egde
                    ds = computeIntegralEdge(i,j,p,q,u,v,U,V,CP,uv,xys);
                    
                    % Insert the load value at the correct location of the
                    % element load vector
                    if direction==1||direction==2
                        Fel(:,:,direction) = h*ds(:,:)*gw*map;
                    elseif direction==3
                        Fel(:,:,1) = h*abs(ds(:,:,2))*gw*map;
                    elseif direction==4
                        Fel(:,:,2) = h*abs(ds(:,:,1))*gw*map;
                    elseif direction==5
                        Fel(:,:,1) = h*ds(:,:,1)*gw*map;
                        Fel(:,:,2) = h*ds(:,:,2)*gw*map;
                    elseif direction==6
                        Fel(:,:,1) = -h*ds(:,:,2)*gw*map;
                        Fel(:,:,2) =  h*ds(:,:,1)*gw*map;
                    end
                    
                    % Assemble the element load vector to the global load
                    % vector in one loop
                    F(i-p:i,j-q:j,:)=Fel(:,:,:)+F(i-p:i,j-q:j,:);
                end  
            end
        end
    end
end 

% Reconstruct the load vector with respect to the global numbering
Fl = make_fl_dof(F);

% If the given old load vector is not null, sum them up with the newly
% computed one
if isvector(Fold)  
    Fl = Fl + Fold;   
end

end