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
