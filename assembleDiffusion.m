function A = assembleDiffusion(mesh, feMap, diffusivity)
    % Gradients of the shape functions in reference coordinates
    shapeGradients = [-1 1 0; -1 0 1];

    % Initialize global stiffness matrix in vector form
    AVector = zeros(9, mesh.numMeshElements);

    % Node indices for assembly
    nodeIndI = [1 2 3 1 2 3 1 2 3];
    nodeIndJ = [1 1 1 2 2 2 3 3 3];

    % Global row and column indices for sparse matrix assembly
    globRows = mesh.meshElements(nodeIndI , :);
    globCols = mesh.meshElements(nodeIndJ , :);

    for e = 1:mesh.numMeshElements
        % Metric tensor for the current element
        C = feMap.metricTensor(:, :, e);
        
        % Local stiffness matrix for the current element with diffusivity
        A_loc = diffusivity(e) * (shapeGradients' * C * shapeGradients / 2);

        % Store the local stiffness matrix in vector form
        AVector(:, e) = A_loc(:);
    end

    % Assemble the global stiffness matrix using sparse format
    A = sparse(globRows(:), globCols(:), AVector(:), mesh.numVertices, mesh.numVertices);
end
