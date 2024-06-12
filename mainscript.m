function runFEMSimulation(meshFile, dt)
    % Load mesh
    mesh = Mesh2D(meshFile);

    % Finite Element Mapping
    feMap = FEMap(mesh);

    % Parameters
    Tf = 35;
    a = 18.515;
    ft = 0.2383;
    fr = 0;
    fd = 1;
    Sigma_h = 9.5298e-4;

    % Initial condition
    u0 = zeros(mesh.numVertices, 1);
    initialNodes = (mesh.vertices(1, :) <= 0.1) & (mesh.vertices(2, :) >= 0.45) & (mesh.vertices(2, :) <= 0.55);
    u0(initialNodes) = 1;

    % Time-stepping
    u = u0;
    numSteps = ceil(Tf / dt);

    % Diffusivity for each element
    Sigma_d = [10 * Sigma_h, Sigma_h, 0.1 * Sigma_h, Sigma_h]; % for different regions
    diffusivity = Sigma_h * ones(mesh.numMeshElements, 1);
    diffusivity(mesh.meshElementFlags == 0) = Sigma_d(1);
    diffusivity(mesh.meshElementFlags == 1) = Sigma_d(2);
    diffusivity(mesh.meshElementFlags == 2) = Sigma_d(3);
    diffusivity(mesh.meshElementFlags == 3) = Sigma_d(4);

    % Assemble mass and stiffness matrices
    M = assembleMass(mesh, feMap);
    K = assembleDiffusion(mesh, feMap, diffusivity);

    % Time integration loop
    activationTimes = zeros(mesh.numVertices, 1);
    activationTimes(:) = Inf; % Set initial activation times to infinity

    for n = 1:numSteps
        % Reaction term
        f = a * (u - fr) .* (u - ft) .* (u - fd);

        % Right-hand side
        rhs = M * u - dt * f;

        % Solve the linear system
        A = M + dt * K;
        u = A \ rhs;

        % Check activation time
        newlyActivated = (u > ft) & (activationTimes == Inf);
        activationTimes(newlyActivated) = n * dt;
    end

    % Plot the final solution
    mesh.plotSolution(u);

    % Check if the matrix is an M-matrix
    isMMatrix = all(all(A <= 0 | (diag(A) > 0)));

    % Check if the potential remains between 0 and 1
    potentialValid = all(u >= 0 & u <= 1);

    % Display results
    fprintf('Results for mesh %s with dt = %f\n', meshFile, dt);
    fprintf('Activation times:\n');
    disp(activationTimes');
    fprintf('Is M-matrix: %d\n', isMMatrix);
    fprintf('Potential within [0, 1]: %d\n', potentialValid);
end

% Function to assemble mass matrix
function M = assembleMass(mesh, feMap)
    % Initialize the mass matrix in vector form
    MVector = zeros(9, mesh.numMeshElements);

    % Shape function values at quadrature points (barycenter in this case)
    shapeValues = [1/3; 1/3; 1/3];
    
    % Jacobian determinant for each element
    J = feMap.J;

    % Node indices for assembly
    nodeIndI = [1 2 3 1 2 3 1 2 3];
    nodeIndJ = [1 1 1 2 2 2 3 3 3];

    % Global row and column indices for sparse matrix assembly
    globRows = mesh.meshElements(nodeIndI, :);
    globCols = mesh.meshElements(nodeIndJ, :);
    
    for e = 1:mesh.numMeshElements
        % Local mass matrix for the current element
        M_loc = (shapeValues * shapeValues') * J(e) / 2;

        % Store the local mass matrix in vector form
        MVector(:, e) = M_loc(:);
    end
    
    % Assemble the global mass matrix using sparse format
    M = sparse(globRows(:), globCols(:), MVector(:), mesh.numVertices, mesh.numVertices);
end

% Function to assemble diffusion matrix
function A = assembleDiffusion(mesh, feMap, diffusivity)
    % Gradients of the shape functions in reference coordinates
    shapeGradients = [-1 1 0; -1 0 1];
    
    % Initialize global stiffness matrix in vector form
    AVector = zeros(9, mesh.numMeshElements);
    
    % Node indices for assembly
    nodeIndI = [1 2 3 1 2 3 1 2 3];
    nodeIndJ = [1 1 1 2 2 2 3 3 3];
    
    % Global row and column indices for sparse matrix assembly
    globRows = mesh.meshElements(nodeIndI, :);
    globCols = mesh.meshElements(nodeIndJ, :);
    
    for e = 1:mesh.numMeshElements
        % Metric tensor for the current element
        C = feMap.metricTensor(:, :, e);
        
        % Local stiffness matrix for the current element
        A_loc = shapeGradients' * C * shapeGradients / 2 * diffusivity(e);
        
        % Store the local stiffness matrix in vector form
        AVector(:, e) = A_loc(:);
    end
    
    % Assemble the global stiffness matrix using sparse format
    A = sparse(globRows(:), globCols(:), AVector(:), mesh.numVertices, mesh.numVertices);
end

% Function to read mesh file
function mesh = readMesh_msh(file_name)
    commandString = ['head -n1 ', file_name, ' > temp.txt'];
    system(commandString);
    data = load('temp.txt');
    
    mesh.number_of_vertices = data(1);
    mesh.number_of_elements = data(2);
    mesh.number_of_boundary_sides = data(3);

    commandString = ['head -n', num2str(data(1) + 1), ' ', file_name, ' > temp.txt'];
    system(commandString);
    data = load('temp.txt');

    mesh.vertices = data(2:end, 1:2)';
    mesh.vertices_flag = data(2:end, 3);

    commandString = ['tail -n', num2str(mesh.number_of_boundary_sides), ' ', file_name, ' > temp.txt'];
    system(commandString);
    data = load('temp.txt');

    mesh.boundary_sides = data(:, 1:2)';
    mesh.boundary_sides_flag = data(:, 3);

    commandString = ['head -n', num2str(mesh.number_of_vertices + mesh.number_of_elements + 1), ' ', file_name, ' > temp.txt'];
    system(commandString);
    commandString = ['tail -n', num2str(mesh.number_of_elements), ' ', 'temp.txt', ' > temp2.txt'];
    system(commandString);
    data = load('temp2.txt');

    mesh.elements = data(:, 1:3)';
    mesh.elements_flag = data(:, 4);

    !rm temp.txt
    !rm temp2.txt

    return
end

% Mesh2D class definition
classdef Mesh2D < handle
    properties
        % Number of vertices
        numVertices
        % Number of mesh elements
        numMeshElements
        % Number of boundary elements
        numBoundaryElements
        % A 2 by numVertices matrix containing coordinates of the vertices
        vertices
        % An array to store flags for vertices
        vertexFlags
        % A 3 by numMeshElements matrix for mesh elements
        meshElements
        % An array to store flags for mesh elements
        meshElementFlags
        % A 2 by numBoundaryElements matrix for boundary elements
        boundaryElements
        % An array to store flags for boundary elements
        boundaryElementFlags
    end
    methods
        % Constructor method to read mesh data from a file
        function obj = Mesh2D(fileName)
            tic
            temp = readMesh_msh(fileName);
            toc

            obj.numVertices = temp.number_of_vertices;
            obj.numMeshElements = temp.number_of_elements;
            obj.numBoundaryElements = temp.number_of_boundary_sides;
            obj.vertices = temp.vertices;
            obj.vertexFlags = temp.vertices_flag;
            obj.meshElements = temp.elements;
            obj.meshElementFlags = temp.elements_flag;
            obj.boundaryElements = temp.boundary_sides;
            obj.boundaryElementFlags = temp.boundary_sides_flag;
        end

        % Method to plot the mesh using trimesh
        function plotMesh(obj)
            close all
            trimesh(obj.meshElements', obj.vertices(1, :), obj.vertices(2, :));
            axis equal
        end

        % Method to plot the mesh with solution values using trimesh
        function plotSolution(obj, solution)
            trimesh(obj.meshElements', obj.vertices(1, :), obj.vertices(2, :), solution);
        end
    end
end

% FEMap class definition
classdef FEMap < handle
    properties
        % The Jacobian of the transformation
        F
        % Translation vector
        b
        % Determinant of the Jacobian
        J % detF
        % Transpose of the Jacobian multiplied by the Jacobian
        C
        % Metric tensor, calculated as J * inv(C)
        metricTensor
    end
    methods
        function obj = FEMap(mesh2D)
            numElements = mesh2D.numMeshElements;

            obj.F = zeros(2, 2, numElements);
            obj.C = zeros(2, 2, numElements);
            obj.metricTensor = zeros(2, 2, numElements);
            obj.b = zeros(2, numElements);
            obj.J = zeros(numElements, 1);

            p1 = mesh2D.vertices(:, mesh2D.meshElements(1, :));
            p2 = mesh2D.vertices(:, mesh2D.meshElements(2, :));
            p3 = mesh2D.vertices(:, mesh2D.meshElements(3, :));

            obj.F(:, 1, :) = p2 - p1;
            obj.F(:, 2, :) = p3 - p1;
            obj.b = p1;

            a = mesh2D.vertices(1, mesh2D.meshElements(1, :));
            b = mesh2D.vertices(1, mesh2D.meshElements(2, :));
            c = mesh2D.vertices(1, mesh2D.meshElements(3, :));
            d = mesh2D.vertices(2, mesh2D.meshElements(1, :));
            e = mesh2D.vertices(2, mesh2D.meshElements(2, :));
            f = mesh2D.vertices(2, mesh2D.meshElements(3, :));

            obj.J = abs(a .* e + b .* f + c .* d - a .* f - b .* d - c .* e);

            for k = 1:numElements
                obj.C(:, :, k) = obj.F(:, :, k)' * obj.F(:, :, k);
                obj.metricTensor(:, :, k) = obj.J(k) * (obj.C(:, :, k) \ eye(2));
            end
        end
    end
end

% Run the simulation for the specified cases
meshes = {'mesh0128.msh', 'mesh0256.msh'};
timeSteps = [0.1, 0.025];

for i = 1:length(meshes)
    for j = 1:length(timeSteps)
        runFEMSimulation(meshes{i}, timeSteps{j});
    end
end
