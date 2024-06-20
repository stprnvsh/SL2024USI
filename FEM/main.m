% List of mesh files and time steps
close all;
clear all; clc;
meshes = {'mesh_0256.msh', 'mesh_0128.msh'}%{'mesh_0064.msh','mesh_0128.msh','mesh_0256.msh'}%, 'mesh_0256.msh','mesh_0064.msh'};
timeSteps = [0.1, 0.05]%, 0.05, 0.025]%, 0.025];
mdiag_flag = 1

% Run simulations for each combination of mesh and time step
for i = 1:length(meshes)
    for j = 1:length(timeSteps)
        runFEMSimulation(meshes{i}, timeSteps(j), mdiag_flag);
    end
end
