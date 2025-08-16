function [dl,bt,sf,ig,names] = define_geometry()

% --- Geometry parameters in a struct ---
ig.numU          = 5;               % number of U-shaped sheets
ig.Nrows         = 4;               % heater rows per U on each inside face
ig.pack_length   = 1;               % depth into page (m)
ig.cell_diameter = 22e-3;           % cell diameter (m)
ig.cell_length   = 70.15e-3;        % cell length (m)
ig.n_cell_module = 18*5;

% U geometry (meters)
ig.a             = 78e-3/2;         % inner half-width of U opening (m)
ig.wall_thk      = 3e-3;             % wall thickness of U (m)
ig.bottom_thk    = 2e-3;             % bottom thickness of U (m)
ig.wall_h        = 120e-3;           % wall height (m)

% baseplate
ig.base_thk      = 4e-3;             % baseplate thickness (m)
ig.margin        = 1e-3;             % side margin around U array (m)

% spacing (center-to-center)
ig.pitch         = 2*(ig.a + ig.wall_thk) + 5e-3;  % adjust gap between Us here

% heater/contact geometry
ig.contact_thk   = 3e-4;             % heater/contact layer thickness (m)
ig.heater_height = ig.cell_diameter; % heater rectangle height

% fin parameters 
ig.Nfins         = 60;               % number of fins
ig.fin_thk       = 2e-3;             % thickness [m]
ig.fin_h         = 5e-3;             % height [m]

% --- Compute U centers along x ---
ig.x_centers = ((0:ig.numU-1) - (ig.numU-1)/2) * ig.pitch;

% prepare gd (10xN rectangle columns), names cell
gd    = zeros(10,0);
names = {};
idx   = 0;

% --- Build U pieces: left wall, right wall, bottom for each U ---
for ui = 1:ig.numU
    xc = ig.x_centers(ui);

    % left wall
    x1 = xc - ig.a - ig.wall_thk; x2 = xc - ig.a;
    y1 = ig.bottom_thk; y2 = ig.wall_h;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_wall_left', ui);

    % right wall
    x1 = xc + ig.a; x2 = xc + ig.a + ig.wall_thk;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_wall_right', ui);

    % bottom (bridge)
    x1 = xc - ig.a - ig.wall_thk; x2 = xc + ig.a + ig.wall_thk;
    y1 = 0; y2 = ig.bottom_thk;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_base', ui);
end

% --- Baseplate rectangle ---
ig.xmin = min(ig.x_centers) - (ig.a + ig.wall_thk) - ig.margin;
ig.xmax = max(ig.x_centers) + (ig.a + ig.wall_thk) + ig.margin;
R = [3;4; ig.xmin; ig.xmax; ig.xmax; ig.xmin; -ig.base_thk; -ig.base_thk; 0; 0];
gd(:,end+1) = R; idx = idx+1; 
names{end+1} = 'baseplate';

% --- Heater/contact strips ---
for ui = 1:ig.numU
    xc = ig.x_centers(ui);
    yc_centers = (1:ig.Nrows) .* (ig.wall_h/(ig.Nrows+1));
    row_h = ig.heater_height;

    for r = 1:ig.Nrows
        yc = yc_centers(r);
        y1 = yc - row_h/2; y2 = yc + row_h/2;

        % left inner face heater
        x1 = xc - ig.a; x2 = xc - ig.a + ig.contact_thk;
        R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
        gd(:,end+1) = R; idx = idx+1; 
        names{end+1} = sprintf('U%d_R%d_thermal_left', ui,r);

        % right inner face heater
        x1 = xc + ig.a - ig.contact_thk; x2 = xc + ig.a;
        R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
        gd(:,end+1) = R; idx = idx+1;
        names{end+1} = sprintf('U%d_R%d_thermal_right', ui,r);
    end
end

% --- Cells inside U channels ---
if isfield(ig,'n_cell_module') && ig.n_cell_module > 0
    cell_w = ig.cell_length;  % width in x-direction
    cell_h = ig.cell_diameter;  % height in y-direction
    
    for ui = 1:ig.numU
        xc = ig.x_centers(ui);  % center x of this U
        yc_centers = (1:ig.Nrows) .* (ig.wall_h/(ig.Nrows+1));  % same vertical spacing as heaters

        for r = 1:ig.Nrows
            yc = yc_centers(r);   % center y position of this row
            
            % rectangle spanning from xc-a → xc+a (inner U walls)
            x1 = xc - ig.a; 
            x2 = xc + ig.a;
            y1 = yc - cell_h/2;
            y2 = yc + cell_h/2;
            
            R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
            gd(:,end+1) = R;
            idx = idx+1;
            names{end+1} = sprintf('U%d_R%d_cell',ui,r);
        end
    end
end

% --- Base plate fins ---
if ig.Nfins > 0
    ig.fin_spacing = (ig.xmax - ig.xmin - ig.Nfins * ig.fin_thk) / (ig.Nfins + 1);
    
    for fi = 1:ig.Nfins
        x1 = ig.xmin + fi * ig.fin_spacing + (fi-1) * ig.fin_thk;
        x2 = x1 + ig.fin_thk;
        y1 = -ig.base_thk - ig.fin_h;
        y2 = -ig.base_thk;
    
        R = [3; 4; x1; x2; x2; x1; y1; y1; y2; y2];
        gd(:, end+1) = R;
        idx = idx + 1;
        names{end+1} = sprintf('fin_%d', fi);
    end
end

% --- Build decsg name arrays & set formula ---
ns = char(names)';           
sf = names{1};
for k = 2:numel(names)
    sf = [sf '+' names{k}];
end

[dl,bt] = decsg(gd,sf,ns);

end