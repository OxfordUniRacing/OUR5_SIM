function [al_faces, ct_faces, cell_faces, cu_faces] = findMaterialFaces(model, ig)
    % FINDMATERIALFACES  Find face IDs for aluminium, TIM, and cells
    %
    % Inputs:
    %   model - PDE ThermalModel
    %   ig    - struct with geometry info (expects .al_centroid, .tim_centroid, .cell_centroid)
    %
    % Outputs:
    %   al_faces   - Face IDs corresponding to aluminium
    %   ct_faces   - Face IDs corresponding to TIM
    %   cell_faces - Face IDs corresponding to cells

    % Aluminium faces
    al_faces = nearestFace(model.Geometry, ig.al_centroid);

    % TIM faces
    ct_faces = nearestFace(model.Geometry, ig.tim_centroid);

    % Cell faces (vector if multiple)
    cell_faces = zeros(size(ig.cell_centroid,1),1);
    for i = 1:size(ig.cell_centroid,1)
        cell_faces(i) = nearestFace(model.Geometry, ig.cell_centroid(i,:));
    end

    % copper faces
    cu_faces = nearestFace(model.Geometry, ig.cu_centroid);
end
