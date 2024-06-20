function M = assembleMass(mesh, feMap)
    % Initialize the mass matrix in vector form
    MVector = zeros(9, mesh.numMeshElements);

    % Shape function values at quadrature points (barycenter in this case)
    shapeValues = [1/3; 1/3; 1/3];
    
    % Jacobian determinant for each element (area * 2 for triangular elements)
    J = feMap.J;

    % Node indices for assembly
    % These indices correspond to the node pairs for each element
    nodeIndI = [1 2 3 1 2 3 1 2 3];
    nodeIndJ = [1 1 1 2 2 2 3 3 3];

    % Global row and column indices for sparse matrix assembly
    % mesh.meshElements(i, :) gives the global node index for the ith local node
    globRows = mesh.meshElements(nodeIndI, :);
    globCols = mesh.meshElements(nodeIndJ, :);
    
    % Loop over each element to compute and assemble local mass matrices
    for e = 1:mesh.numMeshElements
        % Local mass matrix for the current element
        % shapeValues * shapeValues' gives a 3x3 matrix of 1/9s
        % J(e) / 2 gives the area of the element
        M_loc = (shapeValues * shapeValues') * J(e) / 2;

        % Store the local mass matrix in vector form for assembly
        MVector(:, e) = M_loc(:);
    end
    
    % Assemble the global mass matrix using sparse format
    % globRows and globCols provide the global indices for the non-zero entries
    M = sparse(globRows(:), globCols(:), MVector(:), mesh.numVertices, mesh.numVertices);
end
