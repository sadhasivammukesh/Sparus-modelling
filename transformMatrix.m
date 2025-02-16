function Mb = transformMatrix(r, M)
    % Mb: mass matrix expressed in the AUVs buoyancy center,
    % r: vector of coordinates from CO to the buoyancy center of the other part,
    % M: mass matrix in the center of the part.

    H = [eye(3), (S_(r))'; zeros(3, 3), eye(3)];

    Mb = H' * M * H;
end