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
    subplot(1, 2, 1);
    mesh.plotSolution(u0);
    title('Initial Condition');

    % Time-stepping
    u = u0;
    numSteps = ceil(Tf / dt);
    
    figure;
    subplot(1, 3, 1);
    mesh.plotSolution(u0);
    title('Initial Condition');

    % Diffusivity for each element
    Sigma_d = [10 * Sigma_h, Sigma_h, 0.1 * Sigma_h, Sigma_h]; % for different regions
    diffusivity = Sigma_h * ones(mesh.numMeshElements, 1);
    diffusivity(mesh.meshElementFlags == 0) = Sigma_d(1); % Omega_d1
    diffusivity(mesh.meshElementFlags == 1) = Sigma_d(2); % Omega_d2
    diffusivity(mesh.meshElementFlags == 2) = Sigma_d(3); % Omega_d3
    diffusivity(mesh.meshElementFlags == 3) = Sigma_d(4); % Omega_h (healthy tissue)

    % Assemble stiffness matrix
    K = assembleDiffusion(mesh, feMap, diffusivity); % Assembled once, used in A

    % Assemble mass matrix
    M = assembleMass(mesh, feMap); % Assembled once, used in A
    %M = diag(sum(M,2))
    % Precompute the IMEX scheme matrix
    A = M + K*dt; % Precomputed matrix used in time-stepping loop
    
    % Time integration loop
    activationTimes = Inf(mesh.numVertices, 1); % Set initial activation times to infinity
    
    % Initialize data storage
    data = struct();
    data.time = (0:numSteps) * dt;
    data.u = zeros(mesh.numVertices, numSteps + 1);
    data.u(:, 1) = u0;
    data.vertices = mesh.vertices;
    
   % Create a VideoWriter object for MP4 with dynamic filename
    [~, meshName, ~] = fileparts(meshFile); % Extract the mesh name from the file path
    % videoFileName = sprintf('fem_simulation_%s_dt%f_without_diag.mp4', meshName, dt);
    uDataFileName = sprintf('fem_simulation_%s_dt%f_simulation_data_without_diag.mat', meshName, dt);
    % v = VideoWriter(videoFileName, 'MPEG-4');
    % open(v);
    
    for n = 1:numSteps
        % Nonlinear reaction term (explicit part)
        f = a * (u - fr) .* (u - ft) .* (u - fd); % Updated each step
        
        % Right-hand side
        rhs = M * u -  f * dt; % Computed using updated u and f
        
        % Solve the linear system (implicit part)
        u_tnew = A \ rhs; % Solve using precomputed A
        
        % Ensure non-negativity and upper bound of the solution
        u_tnew = max(min(u_tnew, 1), 0);
        
        % Check for NaN values
        if any(isnan(u_tnew))
            fprintf('NaN values detected at step %d\n', n);
            break;
        end

        % Update u
        u = u_tnew; % Updated each step
        
        % Plot intermediate solutions every 100 time steps
        if mod(n, 10) == 0
            mesh.plotSolution(u);
            view(2)
            axis("tight")
            title(['Solution at time step ', num2str(n)]);
            drawnow;
            % Capture the frame for the video
            % frame = getframe(gcf);
            % writeVideo(v, frame);
        end

        % Print intermediate values for debugging
        fprintf('Step %d, max u: %e, min u: %e\n', n, max(u), min(u));

        % Check activation time
        newlyActivated = (u > ft) & (activationTimes == Inf);
        activationTimes(newlyActivated) = n * dt;

        % Store the solution
        data.u(:, n + 1) = u;
    end

    % Plot the final solution
    subplot(1, 2, 2);
    mesh.plotSolution(u);
    title('Final Solution');

    % Check if the potential remains between 0 and 1
    potentialValid = all(u >= 0 & u <= 1);

    % Display results
    fprintf('Results for mesh %s with dt = %f\n', meshFile, dt);
    fprintf('Activation times:\n');
    % disp(activationTimes');
    fprintf('Potential within [0, 1]: %d\n', potentialValid);

    % Save the data to a .mat file
    save(uDataFileName, 'data');
    % figure; % Play the movie at 10 frames per second
    % Close the video file
    % close(v);
end

