function [freeEdges, cellFaces] = getBCindexes(model, ig)

% ============ Cooling edges ============
cooledEdgesCoords = []; 

if ig.Nfins > 0
    % Get all X edges along bottom of baseplate
    y_bottombase = -ig.base_thk;  
    y_midfin     = -ig.base_thk - ig.fin_h/2;
    y_endfin     = -ig.base_thk - ig.fin_h;

    x_start = ig.xmin;
    for fi = 1:ig.Nfins+1
        if fi <= ig.Nfins
            % Gap ends at start of a fin
            x_end = ig.xmin + fi * ig.fin_spacing + (fi-1) * ig.fin_thk;
        else
            x_end = ig.xmax; % After last fin
        end
        
        % Midpoint of this gap segment (bottom edge)
        xm = (x_start + x_end) / 2;
        cooledEdgesCoords(end+1,:) = [xm, y_bottombase];
   
        if fi <= ig.Nfins
            % Midpoint of fin bottom
            x_start_fin = ig.xmin + fi * ig.fin_spacing + (fi-1) * ig.fin_thk;
            x_end_fin   = ig.xmin + fi * ig.fin_spacing + fi * ig.fin_thk;
            xm = (x_end_fin + x_start_fin) / 2;
            cooledEdgesCoords(end+1,:) = [xm, y_endfin];

            % Left and right edges of fin (mid-height)
            cooledEdgesCoords(end+1,:) = [x_start_fin, y_midfin];
            cooledEdgesCoords(end+1,:) = [x_start_fin + ig.fin_thk, y_midfin];

            % Advance start to after this fin
            x_start = x_end + ig.fin_thk;
        end
    end
else
    cooledEdgesCoords = [(ig.xmin + ig.xmax)/2, -ig.base_thk];
end

freeEdges = nearestEdge(model.Geometry, cooledEdgesCoords);


% ============ Cell faces ============
yc_centers = (1:ig.Nrows) .* (ig.module_h / (ig.Nrows + 1));

cellCenters = []; 
for ui = 1:ig.numU
    xc = ig.x_centers(ui);
    for r = 1:ig.Nrows
        yc = yc_centers(r);
        % midpoint of cell region (not heater strip)
        cellCenters(end+1,:) = [xc, yc];
    end
end

cellFaces = nearestFace(model.Geometry, cellCenters);

end