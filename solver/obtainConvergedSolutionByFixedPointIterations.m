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
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function does the fixed point iterations to obtain a solution to
%   the FSI problem.
%
%
%
%
function [str_disp, str_vel, str_acc, ...
          fmeshDisp, fmeshVel, fmeshAcc, ...
          up, upRate, ...
          fFSI, dFSI, ...
          mesh_s_w, mesh_f] = obtainConvergedSolutionByFixedPointIterations(mesh_s,mesh_s_w,bc_s,physics_s,str_disp,str_vel,str_acc, ...
                                                                                        fFSI, fFSIprev, mesh_ale,bc_ale,physics_ale, ...
                                                                                        fmeshDisp, fmeshVel, fmeshAcc, mesh_f,refMesh_f, ...
                                                                                        bc_f,physics_f,up, upRate, fsiProp,transient, t, dFSIVec, fFSIVec)

current_residual     = 10;
iter_count           = 0;
alpha = transient.fsiProp.omega;
% This mesh_s_w is for writing structural result
mesh_s_w            = mesh_s;
% Estimation of initial force
fFSI_fp             = 2.5*fFSIVec(:,3) - 2*fFSIVec(:,2) + 0.5*fFSIVec(:,1);
%%% Estimation of initial displacements and velocities
% Displacements
disp_fsi_prev_fp    = 2.5*dFSIVec(:,3) - 2*dFSIVec(:,2) + 0.5*dFSIVec(:,1);
disp_fsi_prev_time  = disp_fsi_prev_fp;
% Velocities
vel_fsi_prev_fp     = str_vel(1:2*length(mesh_s.nodes));
vel_fsi_prev_fp     = vel_fsi_prev_fp(mesh_s.fsiDOF)';
vel_fsi_prev_fp     = vel_fsi_prev_fp';


%% Starting the Fixed point iteration loop
while(current_residual>fsiProp.fsiResidual && iter_count < fsiProp.maxFpIterations)
    % updating/backup
    up_fp           = up;
    upRate_fp       = upRate;
    
    fprintf('\t \t FIXED POINT ITERATION Number :: %d \n', iter_count);
    %% %%%%%%%%%%%%%%%%% SOLVING STRUCTURAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('######################## SOLVING THE STRUCTURAL PROBLEM #########################\n');
    [disp_str_fp, str_vel_fp, str_acc_fp,dUdot_dU_fsi,petVel] = solve_TransientPlateInMembraneAction(mesh_s, bc_s, physics_s, transient, str_disp, str_vel, str_acc, fFSI_fp);
    % Extracting the velocity and displacements on the FSI interface
    disp_fsi_fp                     = disp_str_fp(mesh_s.fsiDOF);
    vel_fsi_fp                      = str_vel_fp(mesh_s.fsiDOF);
    
    % Creating mesh for writing out structure displacement
    mesh_s_w.nodes(:,1)          = mesh_s_w.nodes(:,1) + disp_str_fp(1:2:end); % X coordinate
    mesh_s_w.nodes(:,2)          = mesh_s_w.nodes(:,2) + disp_str_fp(2:2:end); % Y coordinate
    
    
    %%% We perform the optimization update only when the the two solvers
    %%% are drifting apart, this happens only after first iteration.
    if(iter_count >= 1)
        %% Calculating the residual for Aitken Iterations
        residual        = disp_fsi_prev_fp  - disp_fsi_fp;
        residualVel     = vel_fsi_prev_fp   - vel_fsi_fp;
        %% Obtaining Aitken relaxation factor and relaxing the displacement on the FSI interface
        alpha = alpha_prev_fp*(dot(residue_prev_fp,(residue_prev_fp - residual)))/(norm(residue_prev_fp - residual)^2);
        fprintf('\t \t Relaxation Factor :: %d \n \t \t Residual is :: %d \n \n', alpha, current_residual);
        
        disp_fsi_fp     = disp_fsi_prev_fp  - (alpha)*residual;
        vel_fsi_fp      = vel_fsi_prev_fp   - (alpha)*residualVel;
        %%% Relative Residual calculation
        normDisp_prevTime = norm(disp_fsi_prev_time);
        if(normDisp_prevTime == 0.0)
            current_residual = norm(residual)/1;
        else
            current_residual = norm(residual)/normDisp_prevTime;
        end
        
        % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % %         %%% Obtaining the optimized interface shape by using the Discret Adjoint method.
        % % %         %%% This should also be done for the coordinates and the velocity.
        % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % %         %%% Calculating the objective function
        % % %
        % % %         %%% Obtain the FSI update by node based shape optimization
        % % % 		  %%% EMPIRE communication happens inside this function
        % % %         [cord_fsi_update] = obtainFSIUpdateByFilteringAndOptimizing(dI_dX);
        % % %
        % % %         % Applying the update to displacements.
        % % %         disp_fsi_fp(2:2:end) = disp_fsi_fp(2:2:end) + cord_fsi_update(2:2:end);
    else
        disp_fsi_fp = alpha*disp_fsi_prev_fp    + (1-alpha)*disp_fsi_fp;
        vel_fsi_fp  = alpha*vel_fsi_prev_fp     + (1-alpha)*vel_fsi_fp;
        residual        = disp_fsi_prev_fp  - disp_fsi_fp;        
    end
    
    %% %%%%%%%%%%%%%%%%% SOLVING ALE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adding the interface displacements to the Drichlet boundary
    %%% conditions of the mesh motion solver
    fprintf('######################## MOVING THE MESH #########################\n');
    bc_ale.fsiDisp              = disp_fsi_fp;
    bc_ale.fsiVel               = vel_fsi_fp;
    
    %%% TODO :: Check if the mesh motion solver uses the velocity vel_fsi_fp as boundary condition. If it is using then make it not to use
    [fmeshDisp_fp, fmeshVel_fp, fmeshAcc_fp] = solve_MeshMovement(mesh_ale, bc_ale, physics_ale, transient, fmeshDisp, fmeshVel, fmeshAcc,dUdot_dU_fsi,petVel);
    %%% Extracting the velocity and displacement vectors after mesh
    %%% movement and assigning the mesh velocity
    mesh_f.velocity(1:3:end)    = fmeshVel_fp(1:2:end); % X component
    mesh_f.velocity(2:3:end)    = fmeshVel_fp(2:2:end); % Y component
    %%% Updating the fluid mesh
    mesh_f.nodes(:,1)           = refMesh_f.nodes(:,1) + fmeshDisp_fp(1:2:end); % X coordinate
    mesh_f.nodes(:,2)           = refMesh_f.nodes(:,2) + fmeshDisp_fp(2:2:end); % Y coordinate
    
    %%% Enforcing that the velocity of the fluid on FSI interface to be equal to the
    %%% mesh velocity and the velocity of the material
    bc_f.drichletVector(3*mesh_f.fsiNodes-2) = fmeshVel_fp(bc_ale.fsiDOF(1:2:end)); % X Component
    bc_f.drichletVector(3*mesh_f.fsiNodes-1) = fmeshVel_fp(bc_ale.fsiDOF(2:2:end)); % Y Component     
    %% %%%%%%%%%%%%%%%%% SOLVING FLUID PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('######################## SOLVING THE FLUID PROBLEM #########################\n');
    transient_f = transient;
    dt = transient.dt/transient.noFluidItrFSI_inFPI;
    transient_f.dt = dt;
    for i_f = 1:transient.noFluidItrFSI_inFPI
        [up_fp, upRate_fp, fFSI_fp] = solve_TransientNavierStokesProblemBossakScheme2D(mesh_f, bc_f, physics_f, transient_f, up_fp, upRate_fp);
    end
    
    %%% Scalling the force for making the simulation stable.
    norm_force_before_scalling = norm(fFSI_fp, inf);
    if(t<0.50)
        fFSI_fp = fFSI_fp*(1-cos(pi*t));
    end
    norm_force_after_scalling = norm(fFSI_fp, inf);
    fprintf('\t \t \n Force on FSI interface before scalling :: %d and after scalling :: %d \n', norm_force_before_scalling, norm_force_after_scalling);
    fprintf('\t \t Residual is :: %d \n \n',current_residual);
    
    % Writing the result to vtk files
    %% Writing Fluid Result
    mesh_f.elemOrder = 3;
    FluidTitle = ['Fluid_Result_at_Time_',num2str(transient.i),'_',num2str(iter_count),'.vtk'];
    writeResultToVtk(mesh_f, up_fp, FluidTitle);
    %%% Writing Structure Result
    mesh_s_w.elemOrder = 3;
    StructureTitle = ['Structure_Result_at_Time_',num2str(transient.i),'_', num2str(iter_count),'.vtk'];
    formulatedStrResult = zeros(3*length(mesh_s_w.nodes),1);
    formulatedStrResult(3*mesh_s_w.fsiNodes-2) = fFSI_fp(1:2:end);
    formulatedStrResult(3*mesh_s_w.fsiNodes-1) = fFSI_fp(2:2:end);
    writeResultToVtk(mesh_s_w, formulatedStrResult, StructureTitle);
    clear formulatedStrResult;
    %%% Writing Fluid mesh displacement Result
    mesh_f.elemOrder = 3;
    FluidMeshDispTitle = ['Fluid_Mesh_Disp_Result_at_Time_',num2str(transient.i),'_',num2str(iter_count),'.vtk'];
    formulatedFMeshResult = zeros(3*length(mesh_f.nodes),1);
    formulatedFMeshResult(3*(1:length(mesh_f.nodes))-2) = fmeshDisp_fp(1:2:end);
    formulatedFMeshResult(3*(1:length(mesh_f.nodes))-1) = fmeshDisp_fp(2:2:end);
    writeResultToVtk(mesh_f, formulatedFMeshResult, FluidMeshDispTitle);
    clear formulatedFMeshResult
    
    %% Updation for the next fixed point iteration
    disp_fsi_prev_fp    = disp_fsi_fp;
    vel_fsi_prev_fp     = vel_fsi_fp;
    residue_prev_fp      = residual;
    alpha_prev_fp        = alpha;
    iter_count           = iter_count + 1;
    
    % End of Fixed point Iterations
end
fprintf('\t \t Fixed Point Iterations at Time %d Converged after :: %d Iterations \n \t \t Residual is :: %d \n \n',t, iter_count, current_residual);

% After convergence we assign the result to the output of fixed point
% iterations.
str_disp    = disp_str_fp;
str_vel     = str_vel_fp;
str_acc     = str_acc_fp;

fmeshDisp   = fmeshDisp_fp;
fmeshVel    = fmeshVel_fp;
fmeshAcc    = fmeshAcc_fp;

up      = up_fp;
upRate  = upRate_fp;

fFSI = fFSI_fp;
dFSI = disp_fsi_fp;

% End of the function
end
