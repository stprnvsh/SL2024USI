% List of mesh files and time steps
meshes = {'mesh_0064.msh','mesh_0128.msh', 'mesh_0256.msh'};
timeSteps = [0.1, 0.025];

% Run simulations for each combination of mesh and time step
for i = 1:length(meshes)
    for j = 1:length(timeSteps)
        runFEMSimulation(meshes{i}, timeSteps(j));
    end
end
