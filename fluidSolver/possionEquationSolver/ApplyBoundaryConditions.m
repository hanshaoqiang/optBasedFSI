function [A, b, freeDOFs, drichletDOFs, dValueVec] = ApplyBoundaryConditions(Global_Stiffness_Matrix, Load_Matrix, grid, bc)
%
% This function applies boundary conditions to the Load_Matrix and returns
% the modified Global_Stiffness_Matrix as A and Modified Load_Matrix as b
%
% Boundary Conditions are speicified only in this function
%   1 is Drichlet Boundary Conditions
%   2 is Neumann Boundary Conditons
%
%   left, right, top, bottom indicates the sides of the rectangulr domain.
%
%

b = zeros(grid.numDOF,1);
drichletDOFs = [];
dValueVec = [];
Neumann_values = zeros(grid.numDOF,1);

%% Treating for LEFT boundary

% DRICHLET
if(bc.left.nDrichlet ~= 0)
    for i = 1:bc.left.nDrichlet
        % Check the correctness
        Dnodes = find( grid.nodes(:,1) == grid.xMin & grid.nodes(:,2) >= bc.left.drichlet(i).interval(1) ...
                        & grid.nodes(:,2) <= bc.left.drichlet(i).interval(2) );
        b(2*Dnodes) = bc.left.drichlet(i).value.x;
        b(2*Dnodes-1) = bc.left.drichlet(i).value.y;
        dValueVec = [dValueVec; b(2*Dnodes); b(2*Dnodes-1)];
        drichletDOFs = [drichletDOFs; 2*Dnodes; 2*Dnodes-1];
    end
end

% Neumann
if(bc.left.nNeumann ~= 0)
    for i = 1:bc.left.nNeumann
        % Check the correctness
        Nnodes = find( grid.nodes(:,1) == grid.xMin & grid.nodes(:,2) >= bc.left.neumann(i).interval(1) ...
                        & grid.nodes(:,2) <= bc.left.neumann(i).interval(2) );
        Neumann_values(2*Nnodes) = bc.left.neumann(i).value.x;
        Neumann_values(2*Nnodes-1) = bc.left.neumann(i).value.y;
        
    end
end

%% Treating for RIGHT boundary

% DRICHLET
if(bc.right.nDrichlet ~= 0)
    for i = 1:bc.right.nDrichlet
        % Check the correctness
        Dnodes = find( grid.nodes(:,1) == grid.xMax & grid.nodes(:,2) >= bc.right.drichlet(i).interval(1) ...
            & grid.nodes(:,2) <= bc.right.drichlet(i).interval(2) );
        b(Dnodes) = bc.right.drichlet(i).value.x;
        b(2*Dnodes-1) = bc.right.drichlet(i).value.y;
        dValueVec = [dValueVec; b(2*Dnodes); b(2*Dnodes-1)];
        drichletDOFs = [drichletDOFs; 2*Dnodes; 2*Dnodes-1];
    end
end

% Neumann
if(bc.right.nNeumann ~= 0)
    for i = 1:bc.right.nNeumann
        % Check the correctness
        Nnodes = find( grid.nodes(:,1) == grid.xMax & grid.nodes(:,2) >= bc.right.neumann(i).interval(1) ...
            & grid.nodes(:,2) <= bc.right.neumann(i).interval(2) );
        Neumann_values(2*Nnodes) = bc.right.neumann(i).value.x;
        Neumann_values(2*Nnodes-1) = bc.right.neumann(i).value.y;
    end
end


%% Treating for TOP boundary

% DRICHLET
if(bc.top.nDrichlet ~= 0)
    for i = 1:bc.top.nDrichlet
        % Check the correctness
        Dnodes = find( grid.nodes(:,2) == grid.yMax & grid.nodes(:,1) >= bc.top.drichlet(i).interval(1) ...
                        & grid.nodes(:,1) <= bc.top.drichlet(i).interval(2) );
        b(2*Dnodes) = bc.top.drichlet(i).value.x;
        b(2*Dnodes-1) = bc.top.drichlet(i).value.y;
        dValueVec = [dValueVec; b(2*Dnodes); b(2*Dnodes-1)];
        drichletDOFs = [drichletDOFs; 2*Dnodes; 2*Dnodes-1];
    end
end

% Neumann
if(bc.top.nNeumann~= 0)
    for i = 1:bc.top.nNeumann
        % Check the correctness
        Nnodes = find( grid.nodes(:,2) == grid.yMax & grid.nodes(:,1) >= bc.top.neumann(i).interval(1) ...
                        & grid.nodes(:,1) <= bc.top.neumann(i).interval(2) );
        Neumann_values(2*Nnodes) = bc.top.neumann(i).value.x;
        Neumann_values(2*Nnodes-1) = bc.top.neumann(i).value.y;
    end
end


%% Treating for BOTTOM boundary

% DRICHLET
if(bc.bottom.nDrichlet ~= 0)
    for i = 1:bc.bottom.nDrichlet
        % Check the correctness
        Dnodes = find( grid.nodes(:,2) == grid.yMin & grid.nodes(:,1) >= bc.bottom.drichlet(i).interval(1) ...
                        & grid.nodes(:,1) <= bc.bottom.drichlet(i).interval(2) );
        b(2*Dnodes) = bc.bottom.drichlet(i).value.x;
        b(2*Dnodes-1) = bc.bottom.drichlet(i).value.y;
        dValueVec = [dValueVec; b(2*Dnodes); b(2*Dnodes-1)];
        drichletDOFs = [drichletDOFs; 2*Dnodes; 2*Dnodes-1];
    end
end

% Neumann
if(bc.bottom.nNeumann ~= 0)
    for i = 1:bc.bottom.nNeumann
        % Check the correctness
        Nnodes = find( grid.nodes(:,2) == grid.yMin & grid.nodes(:,1) >= bc.bottom.neumann(i).interval(1) ...
                        & grid.nodes(:,1) <= bc.bottom.neumann(i).interval(2) );
        Neumann_values(2*Nnodes) = bc.bottom.neumann(i).value.x;
        Neumann_values(2*Nnodes-1) = bc.bottom.neumann(i).value.y;
    end
end


%% Treating for the accumulated Drichlet and Neumann values
b = Load_Matrix - Global_Stiffness_Matrix*b;
b = b + Neumann_values;
freeDOFs = setdiff(1:grid.numDOF, drichletDOFs);
b = b(freeDOFs);
A = Global_Stiffness_Matrix(freeDOFs, freeDOFs);
end
