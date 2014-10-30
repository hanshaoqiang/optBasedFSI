function [up] = readInitialConditionsFromVTK(casename)
%% This function reads the output files form the GID output and formulates
% the necessary arrays and element tables for both fluid and structural
% solvers.
% 
% The mesh files are by default assumed to be in the path 
%   INPUT 
%       casename    : name of the vtk file from where the initial values
%                       should be read.
%
%   OUTPUT
%       up          : velocity and pressure at each node.

meshFileName = [casename '.vtk'];
fstring = fileread(meshFileName); % read the mesh file as one string

%% Reading Pressure
fblocks = regexp(fstring,'SCALARS pressure double\nLOOKUP_TABLE default','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
pressure = cell2mat(out(:,1));



%% Reading Velocity
fblocks = regexp(fstring,'VECTORS velocity double','split'); % uses any single character as a separator
fblocks(1) = []; % removes anything before the first character
out = cell(size(fblocks));
for k = 1:numel(fblocks)
    out{k} = textscan(fblocks{k},'%f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
velocity = cell2mat(out(:,1));

% Initializing up vector
up = zeros(3*length(pressure),1);

% Assigning pressure
up(3:3:end) = pressure;

% Assigning velocity
up(1:3:end) = velocity(:,1); % X Component
up(2:3:end) = velocity(:,2); % Y Component


clear pressure velocity



end