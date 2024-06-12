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
