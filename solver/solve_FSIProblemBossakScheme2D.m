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
function solve_FSIProblemBossakScheme2D(mesh_f, bc_f, physics_f, mesh_s, bc_s, physics_s, mesh_ale, bc_ale, physics_ale, transient)
%% Function documentation
%
% Solve the FSI problem in 2D using the Bossak time
% integration scheme for the temporal discretization and the Triangular basis for the spatial discretization.
%
% Unconditionally stability is ensured if the following relations hold :
%
% - alphaBeta <= .5
% - beta >= gamma/2 >= .25
% - alphaBeta + gamma >= .25
%
%   INPUT:
%       mesh_f                  :   mesh of the fluid containing the mesh(nodes, elements, boundaryElements ...)
%       bc_f                    :   boundary conditions applied on the
%                                   fluid formulated in terms of nodes and
%                                   vectors.
%       physics_f              :   physicsl properties of the fluid like
%                                   the kinematic viscosity, nonlinear solver properties.
%
%
%       mesh_s                  :   mesh of the structure containing the mesh(nodes, elements, boundaryElements ...)
%       bc_s                    :   boundary conditions applied on the
%                                   structure formulated in terms of nodes and
%                                   vectors.
%       physics_s               :   physicsl properties of the structure like
%       mesh_ale                :   mesh of the ale fluid domain containing the mesh(nodes, elements, boundaryElements ...)
%       bc_ale                  :   boundary conditions applied on the
%                                   ale fluid mesh formulated in terms of nodes and
%                                   vectors.
%       physics_ale             :   physicsl properties of the ale like
%                                   the Youngs modulus, tolerance etc.
%       transient               :   parameters for time integration scheme
%
%   OUTPUT:
%       upHistory               :   Solution of the Navier-Stokes problem on fluid domain at
%                                   every time step
%       dispHistory             :   Displacement solution of the structural domain at
%                                   every time step
%
%% Function main body
fprintf('__________________________________________________\n');
fprintf('\n');
fprintf('Nonlinear transient analysis for FSI problem in 2D \n');
fprintf('using the Bossak time integration scheme for both fluid and \n');
fprintf('structural domain has been initiated \n');
if strcmp(physics_f.nonLinearScheme.method,'Newton')
    fprintf('Applied nonlinear scheme: Newton method \n');
    fprintf('Residual tolerance: %d \n',physics_f.nonLinearScheme.eps);
    fprintf('Maximum number of iterations: %d \n',physics_f.nonLinearScheme.maxIter);
    fprintf('Residual tolerance for fixed point iterations: %d \n',transient.fsiProp.fsiResidual);
    fprintf('Maximum number of FSI fixed point iterations: %d \n',transient.fsiProp.maxFpIterations);
end
fprintf('__________________________________________________\n');
fprintf('\n');


%% Initializing the initial solution vectors and other required vectors

%%%%%%%%%%%%%%%%%%%%% FLUID %%%%%%%%%%%%%%%%%%%%%%%%
% Compute the number of degrees of freedom
fnDoFs = 3*length(mesh_f.nodes); % number of DOF for fluid solver
F = zeros(fnDoFs,1); % Initialize the global flux vector
%up = zeros(length(F),1);% Get the initial condition for the pressure and velocity field 
up = readInitialConditionsFromVTK('Initial'); % up = ones(length(F),1);
upRate = zeros(length(F),1);% Assume that the pressure and velocity rate is zero initially

%%%%%%%%%%%%%%%%%%%%% FLUID MESH %%%%%%%%%%%%%%%%%%%%
disp_fmesh = zeros(2*length(mesh_f.nodes),1);
vel_fmesh = zeros(2*length(mesh_f.nodes),1);
acc_fmesh = zeros(2*length(mesh_f.nodes),1);
mesh_f.velocity = zeros(3*length(mesh_f.nodes),1);
refMesh_f = mesh_f; % Storing the mesh which will be updated

%%%%%%%%%%%%%%%%%%%%% STRUCTURE %%%%%%%%%%%%%%%%%%%%
snDoFs = 2*length(mesh_s.nodes); % Compute the number of degrees of freedom
% Assuming the initial displacement, velocity and accelaration of the structure to be zero
str_disp = zeros(snDoFs,1); % Displacement
str_vel = zeros(snDoFs,1);  % Velocity
str_acc = zeros(snDoFs,1);  % Acceleration

%%%%%%%%%%%%%%%%%%%%% FSI Related %%%%%%%%%%%%%%%%%%%%%%
force_fsi = zeros(length(mesh_s.fsiDOF),1); % Force on the FSI interface
fFSIprev = zeros(length(mesh_s.fsiDOF),1); % Force on the FSI interface
dFSIVec = zeros(length(mesh_s.fsiDOF),3); % vector of displacements on the FSI interface
fFSIVec = zeros(2*length(mesh_s.fsiNodes),3); 
%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%
% Initialize the simulation time
t = transient.Tstart;

%%%%%%%%%%%%%%%%%%%%% EMPIRE CONNECTION  %%%%%%%%%%%%%%%%%%%%%%


% The time step for the transient analysis
transient.i = 1;
count = 1;
%% Loop over all the time instances of the simulation
fprintf('\t Looping over all the time steps \n');
fprintf('\t ------------------------------- \n');
fprintf('\n');
while(t<=transient.Tend)
    
    mesh_s_w = mesh_s;
    timeTitle = sprintf(' \n \n \n \t \t ::::::: TIME t = %d ::::::: \n \n',t);
    fprintf(timeTitle);
    % Starting fixed point iterations
    if(transient.i > 1)      
        %% Performing the fixed point iterations to obtain the converged solution
        % This method solves for strucuture,
        % fluid, and mesh movement in every fixed point iteration and
        % returns the pressure-velocity and displacement fields of
        % fluid and structure respectively.
        % Mesh is also moved accordingly
        [str_disp, str_vel, str_acc, ...
         disp_fmesh, vel_fmesh, acc_fmesh, ...
         up, upRate, force_fsi, disp_fsi, ...
         mesh_s_w, mesh_f]       = obtainConvergedSolutionByFixedPointIterations(mesh_s, mesh_s_w, bc_s, ...
                                                                                 physics_s, str_disp, str_vel, ...
                                                                                 str_acc, force_fsi, fFSIprev, ...
                                                                                 mesh_ale, bc_ale, physics_ale,...
                                                                                 disp_fmesh, vel_fmesh, acc_fmesh, ...
                                                                                 mesh_f, refMesh_f, bc_f, physics_f, up, ...
                                                                                 upRate, transient.fsiProp, transient, t, dFSIVec, fFSIVec);
                                                  
         
        count = count + 1;        
    end
    
    %% %%%%%%%%%%%%%%%%% SOLVING FLUID PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(transient.i <= 1)
        fprintf('######################## SOLVING THE FLUID PROBLEM #########################\n');
        
        transient_f = transient;
        noFluidItrFSI = transient.noFluidItrFSI;
        
        for i_f = 1:noFluidItrFSI
            [up, upRate, force_fsi] = solve_TransientNavierStokesProblemBossakScheme2D(mesh_f, bc_f, physics_f, transient_f, up, upRate);
        end
        
        % Scalling the force acting on FSI interface
        if(t<0.50)
            force_fsi = force_fsi*(1-cos(pi*t));
        end
        disp_fsi = zeros(length(mesh_s.fsiDOF),1);
    end
    
    
    % Saving the FSI displacements
    dFSIVec(:,1) = dFSIVec(:,2); 
    dFSIVec(:,2) = dFSIVec(:,3);
    dFSIVec(:,3) = disp_fsi;
    
    
    fFSIVec(:,1) = fFSIVec(:,2); 
    fFSIVec(:,2) = fFSIVec(:,3);
    fFSIVec(:,3) = force_fsi; 
    
      
    %% Visualization
    if(mod(transient.i,1) == 0)
        %%% Writing Fluid Result
        mesh_f.elemOrder = 3;
        FluidTitle = ['Fluid_Result_at_Time_',num2str(transient.i),'.vtk'];
        writeResultToVtk(mesh_f, up, FluidTitle);
        %%% Writing Structure Result
        mesh_s_w.elemOrder = 3;
        StructureTitle = ['Structure_Result_at_Time_',num2str(transient.i),'.vtk'];
        formulatedStrResult = zeros(3*length(mesh_s_w.nodes),1);
        formulatedStrResult(3*mesh_s_w.fsiNodes-2) = force_fsi(1:2:end);
        formulatedStrResult(3*mesh_s_w.fsiNodes-1) = force_fsi(2:2:end);
        writeResultToVtk(mesh_s_w, formulatedStrResult, StructureTitle);
        clear formulatedStrResult;
    end
    
    %%% updating the time and other entities for next time step
    t = t + transient.dt;
    transient.i = transient.i + 1;
    
    % End of time loop
end

% End of function
end
