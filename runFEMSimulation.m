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

    % Plot the initial conditions
    figure;
    subplot(1, 3, 1);
    mesh.plotSolution(u0);
    title('Initial Condition');

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

    % Plot the diffusivity setup
    figure;
    hold on;
    colormap(jet(4));
    for i = 1:mesh.numMeshElements
        if mesh.meshElementFlags(i) == 0
            color = 'r'; % Region 1
        elseif mesh.meshElementFlags(i) == 1
            color = 'g'; % Region 2
        elseif mesh.meshElementFlags(i) == 2
            color = 'b'; % Region 3
        else
            color = 'y'; % Healthy tissue
        end
        fill(mesh.vertices(1, mesh.meshElements(:, i)), mesh.vertices(2, mesh.meshElements(:, i)), color);
    end
    hold off;
    title('Diffusivity Regions');
    axis equal;
    colorbar;
    caxis([0, 3]);
    legend('Region 1', 'Region 2', 'Region 3', 'Healthy tissue');

    % Assemble mass and stiffness matrices
    M = assembleMass(mesh, feMap);
    M_lumped = diag(sum(M, 2));
    K = assembleDiffusion(mesh, feMap, diffusivity);

    % Precompute the IMEX scheme matrices
    A = M - dt * K;
    A_lumped = M_lumped - dt * K;

    % Check if A is an M-matrix (efficiently)
    diagA = spdiags(A, 0);
    offDiagA = A - spdiags(diagA, 0, size(A, 1), size(A, 2));
    isMMatrix = all(diagA > 0) && nnz(offDiagA > 0) == 0;

    % Time integration loop
    activationTimes = Inf(mesh.numVertices, 1); % Set initial activation times to infinity

    use_lumped_mass = false; % Flag to switch to lumped mass matrix

    for n = 1:numSteps
        % Reaction term
        f = a * (u - fr) .* (u - ft) .* (u - fd);

        % Right-hand side
        if use_lumped_mass
            rhs = M_lumped * u - dt * f;
        else
            rhs = M * u - dt * f;
        end

        % Solve the linear system
        if use_lumped_mass
            u = A_lumped \ rhs;
        else
            u = A \ rhs;
        end

        % Check activation time
        newlyActivated = (u > ft) & (activationTimes == Inf);
        activationTimes(newlyActivated) = n * dt;

        % Ensure potential remains within bounds
       
        % Check if potential exceeds bounds and switch to lumped mass matrix if needed
        if any(u < 10^-10) || any(u > 1 + 10^-10)
            use_lumped_mass = true;
        end
    end

    % Plot the final solution
    subplot(1, 3, 2);
    mesh.plotSolution(u);
    title('Final Solution');

    % Check if the potential remains between 0 and 1
    potentialValid = all(u >= 0 & u <= 1);

    % Display results
    fprintf('Results for mesh %s with dt = %f\n', meshFile, dt);
    fprintf('Is M-matrix: %d\n', isMMatrix);
    fprintf('Potential within [0, 1]: %d\n', potentialValid);
end
