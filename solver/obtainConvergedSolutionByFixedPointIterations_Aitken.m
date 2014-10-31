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
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sdisp, Svel, Sacc, fmeshDisp, fmeshVel, fmeshAcc, up, upRate, fFSI, dFSI, mesh_s_w, mesh_f] = obtainConvergedSolutionByFixedPointIterations_Aitken(mesh_s,mesh_s_w,bc_s,physics_s,Sdisp,Svel,Sacc,fFSI, fFSIprev, ...
    mesh_ale,bc_ale,physics_ale,fmeshDisp, fmeshVel, fmeshAcc, ...
    mesh_f,refMesh_f,bc_f,physics_f,up, upRate, fsiProp,transient, t, dFSIVec)
%   This function does the fixed point iterations using the steady state
%   solutions using the provided values as initial values

residual = 1;
iterCount = 0;
fsiDisp_prevTime = Sdisp(1:2*length(mesh_s.nodes));
fsiDisp_prevTime = fsiDisp_prevTime(mesh_s.fsiDOF)';
fFSI_fp = fFSI;
fsiDispPrev_fp = 2.5*dFSIVec(:,3) - 2*dFSIVec(:,2) + 0.5*dFSIVec(:,1);
fsiDispPrev_fp = fsiDispPrev_fp';

fsiVel_prevTime = Svel(1:2*length(mesh_s.nodes));
fsiVel_prevTime = fsiVel_prevTime(mesh_s.fsiDOF)';
fsiVel_fp = fsiVel_prevTime';

numFSINodes = length(fFSI)/2;

fFSIPrev_fp = fFSI;
up_fp = up;
upRate_fp = upRate;

% Aitken Relaxation factor
alpha = fsiProp.omega;
normRes = 1;
currentResidual = 1;
initialResidual = 1;
msgEmpty = '';
dt = transient.dt;

%% Starting the Fixed point iteration loop
while(currentResidual>fsiProp.fsiResidual && iterCount < fsiProp.maxFpIterations)
    
    up_fp = up;
    upRate_fp = upRate;
    
    Sdisp_fp = Sdisp;
    Svel_fp = Svel;
    Sacc_fp = Sacc;
    

    mesh_s_w = mesh_s;
    reverseStrFP ='';
    strFpTitle = sprintf(' \t \t FIXED POINT ITERATION Number :: %d \n', iterCount);
    fprintf([reverseStrFP, strFpTitle]);
    
    if(transient.fsiProp.aitkenResidualType == 1 || iterCount > 0)
        
        %% %%%%%%%%%%%%%%%%% SOLVING STRUCTURAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%
        reverseStrStru = '';
        strSolveTitle = sprintf('######################## SOLVING THE STRUCTURAL PROBLEM ######################### \n');
        fprintf([reverseStrStru, strSolveTitle]);
        [Sdisp_fp, Svel_fp, Sacc_fp] = solve_TransientPlateInMembraneAction(mesh_s, bc_s, physics_s, transient, Sdisp, Svel, Sacc, fFSI_fp);
        % Extracting the displacement and velocity from the solution
        disp_s = Sdisp_fp;
        vel_s = Svel_fp;
        % Extracting the velocity and displacements on the FSI interface
        fsiDisp = disp_s(mesh_s.fsiDOF)';
        fsiVel  = vel_s(mesh_s.fsiDOF);
        
        % Creating mesh for writing out structure displacement
        mesh_s_w.nodes(:,1) = mesh_s_w.nodes(:,1) + disp_s(1:2:end); % X coordinate
        mesh_s_w.nodes(:,2) = mesh_s_w.nodes(:,2) + disp_s(2:2:end); % Y coordinate
        
        if(transient.fsiProp.aitkenResidualType == 2)
            %% Calculating the residual for Aitken Iterations
            % Force based residual
            residual = fsiDispPrev_fp - fsiDisp;
            residualVel = fsiVel_fp - fsiVel;
            % Calculating Aitken relaxation factor only from second iteration
            if(iterCount > 1)
                alpha = alphaPrev_fp*(dot(residuePrev_fp,(residuePrev_fp - residual)))/(norm(residuePrev_fp - residual)^2);
            end
            
            %% Obtaining Aitken relaxation factor and relaxing the displacement on the FSI interface
            if(iterCount > 1)
                fsiDisp = fsiDispPrev_fp - (alpha)*residual;
                fsiVel = fsiVel_fp - alpha*residualVel;
            else
                fsiDisp = alpha*fsiDispPrev_fp + (1-alpha)*fsiDisp;
                fsiVel = alpha*fsiVel_fp + (1-alpha)*fsiVel;
            end
            normDisp_prevTime = norm(fsiDisp_prevTime);
            if(normDisp_prevTime == 0.0)
                currentResidual = norm(residual)/1;
            else
                currentResidual = norm(residual)/normDisp_prevTime;
            end
            
            % Updating the iteration count
            fsiDispPrev_fp = fsiDisp;
            fFSIPrev_fp = fFSI_fp;
            fsiVel_fp = fsiVel;
            residuePrev_fp = residual;
            alphaPrev_fp = alpha;
            fprintf('\t \t Relaxation Factor :: %d \n \t \t Residual is :: %d \n \n', alpha, currentResidual);
        end
        
        %% %%%%%%%%%%%%%%%%% SOLVING ALE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adding the interface displacements to the Drichlet boundary
        % conditions of the mesh motion solver
        reverseStrMesh = '';
        strMeshMoveTitle = sprintf('######################## MOVING THE MESH ######################### \n');
        fprintf([reverseStrMesh, strMeshMoveTitle]);
        bc_ale.fsiDisp = fsiDisp;
        bc_ale.fsiVel = fsiVel;
        
        [fmeshDisp_fp, fmeshVel_fp, fmeshAcc_fp] = solve_MeshMovement(mesh_ale, bc_ale, physics_ale, transient, fmeshDisp, fmeshVel, fmeshAcc);
        % Extracting the velocity and displacement vectors after mesh
        % movement
        FmeshDisp =  fmeshDisp_fp;
        FmeshVel =  fmeshVel_fp;
        % Assigning the mesh velocity
        mesh_f.velocity(1:3:end) = FmeshVel(1:2:end); % X component
        mesh_f.velocity(2:3:end) = FmeshVel(2:2:end); % Y component
        % Updating the fluid mesh
        mesh_f.nodes(:,1) = refMesh_f.nodes(:,1) + FmeshDisp(1:2:end); % X coordinate
        mesh_f.nodes(:,2) = refMesh_f.nodes(:,2) + FmeshDisp(2:2:end); % Y coordinate
        
        % Enforcing that the velocity of the fluid on FSI interface to be equal to the
        % mesh velocity and the velocity of the material
        bc_f.drichletVector(3*mesh_f.fsiNodes-2) = fsiVel(1:2:end); % X Component
        bc_f.drichletVector(3*mesh_f.fsiNodes-1) = fsiVel(2:2:end); % Y Component
    
    end
    %% %%%%%%%%%%%%%%%%% SOLVING FLUID PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reverseStrFluid = '';
    strFluidSolverTitle = sprintf('######################## SOLVING THE FLUID PROBLEM #########################\n');
    fprintf([reverseStrFluid, strFluidSolverTitle]);
    transient_f = transient;
    transient_f.dt = dt/transient.noFluidItrFSI_inFPI;
    
    for i_f = 1:transient.noFluidItrFSI_inFPI
        [up_fp, upRate_fp, fFSI_fp] = solve_TransientNavierStokesProblemBossakScheme2D(mesh_f, bc_f, physics_f, transient_f, up_fp, upRate_fp);
    end
    
    fprintf('\t \t \n Inside FIXED PT ITER -- Magnitude of Force on FSI interface BEFORE scaling :: %d \n',norm(fFSI_fp,inf));
    if(t<0.50)
        fFSI_fp = fFSI_fp*(1-cos(pi*t));
    end
    fprintf('\t \t \n Inside FIXED PT ITER -- Magnitude of Force on FSI interface :: %d \n',norm(fFSI_fp,inf));
    
    reverseStrFP = repmat(sprintf('\b'), 1, length(strFpTitle));
    fprintf([reverseStrFP, msgEmpty]);
    
    
    if(transient.fsiProp.aitkenResidualType == 1)
        %% Calculating the residual for Aitken Iterations
        % Force based residual
        residual = fFSIPrev_fp - fFSI_fp;
        % Calculating Aitken relaxation factor only from second iteration
        if(iterCount > 1 && norm(residual,inf))
            alpha = alphaPrev_fp*(dot(residuePrev_fp,(residuePrev_fp - residual)))/(norm(residuePrev_fp - residual)^2);
        end
        
        %% Obtaining Aitken relaxation factor and relaxing the displacement on the FSI interface
        if(iterCount > 1)
            fFSI_fp = fFSIPrev_fp - (alpha)*residual;
        else
            fFSI_fp = alpha*fFSIPrev_fp + (1-alpha)*fFSI_fp;
        end
        
        
        % Updating the iteration count
        fsiDispPrev_fp = fsiDisp;
        fFSIPrev_fp = fFSI_fp;
        fsiVel_fp = fsiVel;
        residuePrev_fp = residual;
        alphaPrev_fp = alpha;
        
        currentResidual = norm(residual)/norm(fFSI);
        fprintf('\t \t Relaxation Factor :: %d \n \t \t Residual is :: %d \n \n', alpha, currentResidual);
       
    end
    
    
    iterCount = iterCount + 1;
    % End of Fixed point Iterations
end

fprintf(' \t \t Fixed Point Iterations at Time %d Converged after :: %d Iterations \n \t \t Residual is :: %d \n \n',t, iterCount, currentResidual);

% After convergence we assign the result to the output of fixed point
% iterations.
Sdisp = Sdisp_fp;
Svel = Svel_fp;
Sacc = Sacc_fp;
fmeshDisp = fmeshDisp_fp;
fmeshVel = fmeshVel_fp;
fmeshAcc = fmeshAcc_fp;
up = up_fp;
upRate = upRate_fp;
fFSI = fFSI_fp;
dFSI = fsiDisp;

% End of the function
end