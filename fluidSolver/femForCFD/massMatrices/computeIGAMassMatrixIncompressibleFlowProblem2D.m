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
function M = computeIGAMassMatrixIncompressibleFlowProblem2D(p,U,q,V,CP,int)
%% Function documentation
%
% Returns the master mass matrix corresponding to an incompressible flow in
% 2D, i.e. there are both velocities and pressure degrees of freedom.
%
%         Input : 
%           p,q : The polynomial degrees of the NURBS patch
%           U,V : The knot vectors of the NURBS patch
%            CP : The Control Points of the NURBS patch
%           int : Structure responsible for the integration
%                     automatic : Automatic choice of the Gauss Points 
%                                 fulfiling stability
%                        manual : Manual choice of the Gauss Points for 
%                                 improving the accuracy
%
%       Output :
%            M : master mass matrix
%
% Function layout :
%
% 0. Read input
%
% 1. Assign DoF numbering
%
% 2. Choose an integration rule
%
% 3. loop over all elements (knot spans)
%
%    3i. Initialize element matrices and force vectors
%
%   3ii. Loop over all Quadrature Points for the integration of the element mass
%
%        3ii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%        3ii.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
%
%        3ii.3. Find the correct spans where u,v lie in
%
%        3ii.4. Compute the NURBS basis functions and their derivatives at the quadrature point
%
%        3ii.5. Compute the determinant of the Jacobian (and possibly Hessian) to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%        3ii.6. Compute the element mass matrix at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
%
% 3iii. Create an element freedom table
%
%  3iv. insert Ke into K via the element freedom table
%
%% Function main body

%% 0. Read input

% Check input
mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
checkBSplineInput2D(p,mu,nu,q,mv,nv);

% Total number of degrees of freedom for the 2D vector transport equation
nDoFs = 3*nu*nv;

% Number of Control Points that affect each knot span
ne = (p+1)*(q+1);

% Number of DoFs affected the element under study
ndof_e = 3*ne;

% Initialize global stiffness matrix
M = zeros(nDoFs,nDoFs);

%% 1. Assign DoF numbering

% dof array assigns two DoF (x,y) to every CP
% numbering follows CP: CP1->dof1,dof2 CP2->dof3,dof4
% Initialize the array of the degrees of freedom
dof = zeros(nu,nv,2);

% Initialize counter
k = 1;

% Loop over all the Control points
for cpj = 1:nv
    for cpi = 1:nu
        dof(cpi,cpj,1) = k;
        dof(cpi,cpj,2) = k + 1;
        dof(cpi,cpj,3) = k + 2;
        
        % Update counter
        k = k + 3;
    end
end

%% 2. Choose an integration rule

% Select the given integration scheme and recover the analogous Gauss
% points :
% default is FGI integration element-wise
if int.type == 0
    nGPu = ceil(p + .5);
    nGPv = ceil(q + .5);
    [GPu,GWu] = gauss(nGPu);
    [GPv,GWv] = gauss(nGPv);
elseif int.type == 1
    [GPu,GWu] = gauss(int.uNGaussMass);
    [GPv,GWv] = gauss(int.vNGaussMass);
end

%% 3. loop over all elements (knot spans)
for j = (q+1):(mv-q-1)   
    for i = (p+1):(mu-p-1)
        % check if we are in a non-zero knot span
        if (U(i+1)~=U(i) && V(j+1)~=V(j))
            %% 3i. Initialize element matrices and force vectors
            
            % Initialize element stiffness matrix
            Me = zeros(ndof_e,ndof_e);
         
            %% 3ii. Loop over all Quadrature Points for the integration of the element mass
            for kv = 1:length(GPv)
                for ku =1:length(GPu)
                    %% 3ii.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    if int.type == 0 || int.type == 1
                        u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
                        v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
                    end
                
                    %% 3ii.2. Issue the quadrature weights the compute the quadrature weight needed in the multiplication
                    gw = GWu(ku)*GWv(kv);
                
                    %% 3ii.3. Find the correct spans where u,v lie in
                    spanu = findKnotSpan(u,U,nu);
                    spanv = findKnotSpan(v,V,nv);
                
                    %% 3ii.4. Compute the NURBS basis functions and their derivatives at the quadrature point
                    [R,dR] = computeNurbsBasisFunctionsAndFirstDerivatives2D(spanu,p,u,U,spanv,q,v,V,CP);
                        
                    %% 3ii.5. Compute the determinant of the Jacobian (and possibly Hessian) to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
                
                    % Initialize Jacobian
                    Jxxi = zeros(2,2);
                
                    % initialize counter
                    k = 0;
                
                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:q
                        for b = 0:p
                            % Update counter
                            k = k + 1;
                        
                            % Compute recursively the entries of the Jacobian
                            Jxxi(1,1) = Jxxi(1,1) + CP(i-p+b,j-q+c,1)*dR(k,1);
                            Jxxi(1,2) = Jxxi(1,2) + CP(i-p+b,j-q+c,2)*dR(k,1);
                            Jxxi(2,1) = Jxxi(2,1) + CP(i-p+b,j-q+c,1)*dR(k,2);
                            Jxxi(2,2) = Jxxi(2,2) + CP(i-p+b,j-q+c,2)*dR(k,2);
                        end
                    end
                
                    % Compute the determinant of the Jacobian
                    detJxxi = det(Jxxi);
                
                    % Compute the determinant of the Jacobian to the
                    % transformation from the NURBS space (xi-eta) to the
                    % integration domain [-1,1]x[-1,1] 
                    %
                    %         | xi_i+1 - xi_i                    |
                    %         | -------------            0       |
                    %         |        2                         |
                    %  xi,u = |                                  |
                    %         |                  eta_j+1 - eta_j |
                    %         |        0         --------------- |
                    %         |                          2       |
                    detJxiu = (U(i+1)-U(i))*(V(j+1)-V(j))/4;
                    
                    %% 3ii.6. Compute the element mass matrix at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
                    
                    % Form the basis function matrix (only zeros will be 
                    % placed in the pressure DoFs since those involve no 
                    % inertia)
                    Rmatrix = zeros(3,ndof_e);

                    % Loop over all the entries of the basis functions matrix
                    for rCounter=1:ne
                        Rmatrix(1,3*rCounter-2) = R(rCounter);
                        Rmatrix(2,3*rCounter-1) = R(rCounter);
                    end
                    
                    MeOnGP = Rmatrix'*Rmatrix;
                    
                    % Add the contribution from the Gauss Point
                    Me = Me + MeOnGP*abs(detJxxi)*abs(detJxiu)*gw;
                end
            end
    
            %% 3iii. Create an element freedom table
            
            % Initialize element freedom table
            dof_l = zeros(1,ndof_e);
            
            % Initialize counter
            k=1;
            
            % Relation global-local DoFs
            for cpj = j-q:j
                for cpi = i-p:i
                    dof_l(k)   = dof(cpi,cpj,1);
                    dof_l(k+1) = dof(cpi,cpj,2);
                    dof_l(k+2) = dof(cpi,cpj,3);
                    
                    % update counter
                    k = k + 3;
                end
            end
            
            %% 3iv. insert Ke into K via the element freedom table
            for kej = 1:3*(1+p)*(1+q)
                % Assemble to the global stiffness matrix
                for kei = 1:3*(1+p)*(1+q)
                    M(dof_l(kei),dof_l(kej)) = M(dof_l(kei),dof_l(kej)) + Me(kei,kej);
                end
            end
        end
    end
end

end
  