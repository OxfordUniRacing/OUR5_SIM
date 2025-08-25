function [dl,bt,sf,ig,names] = define_geometry_2()

% --- Geometry parameters in a struct ---
ig.numU          = 5;               % number of U-shaped sheets
ig.Nrows         = 5;               % heater rows per U on each inside face
ig.pack_length   = 0.5;             % module depth into page (m)
ig.cell_diameter = 22e-3;           % cell diameter (m)
ig.cell_length   = 70.15e-3;        % cell length (m)
ig.n_cell_module = 18*5;

% thermal pad geometry
ig.contact_thk   = 3e-3;             % heater/contact layer thickness (m)
ig.heater_height = ig.cell_diameter; % heater rectangle height

% U geometry (meters)
ig.a             = ig.cell_length / 2;         % inner half-width of U opening (m)
ig.busbar_thk    = 3e-3;             % wall thickness of U (m)
ig.bottom_thk    = ig.busbar_thk;             % bottom thickness of U (m)
ig.module_h        = 120e-3;           % wall height (m)
ig.NbusbarModule = 18; % number of busbars per module
ig.busbar_width = 25e-3;

% baseplate
ig.base_thk      = 4e-3;             % baseplate thickness (m) (4 mm thick base)
ig.margin        = 1e-3;             % side margin around U array (m)

% spacing (center-to-center)
ig.pitch         = 2*(ig.a + ig.busbar_thk) + 5e-3;  % adjust gap between Us here


% fin parameters (assuming using ATS-EXL78-1220-R0 cut down into sections across base, https://www.qats.com/Product/Heat-Sinks/Extrusion-Profiles-lengths-/Profiles/ATS-EXL78-1220-R0/3664.aspx)

ig.Nfins         = 0;% floor(24/146e-3 * 400e-3);% number of fins (heatsink has 24 fins per 146mm width, cut down to fit over approximate base width)
ig.fin_thk       = 0.5e-3;             % thickness [m] (approxiamate)
ig.fin_h         = 4.6e-3;             % height [m]
if ig.Nfins>0
     ig.bottom_thk = ig.bottom_thk + 3.8e-3; % plus 3.8mm for heatsink base
end

% --- Compute U centers along x ---
ig.x_centers = ((0:ig.numU-1) - (ig.numU-1)/2) * ig.pitch;

% prepare gd (10xN rectangle columns), names cell
gd    = zeros(10,0);
names = {};
idx   = 0;
ig.cu_centroid   = [];   % aluminium walls + base
ig.al_centroid   = [];   % aluminium walls + base
ig.tim_centroid  = [];   % thermal pads
ig.cell_centroid = [];   % cells

% --- Build U pieces: left wall, right wall, bottom for each U ---
for ui = 1:ig.numU
    xc = ig.x_centers(ui);

    % left busbar
    x1 = xc - ig.a - ig.busbar_thk; x2 = xc - ig.a;
    y1 = ig.bottom_thk + ig.contact_thk; y2 = ig.module_h + ig.contact_thk;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_busbar_left', ui);
    ig.cu_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];

    % right busbar
    x1 = xc + ig.a; x2 = xc + ig.a + ig.busbar_thk;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_busbar_right', ui);
    ig.cu_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];

    % bottom (busbar)
    x1 = xc - ig.a - ig.busbar_thk; x2 = xc + ig.a + ig.busbar_thk;
    y1 = ig.contact_thk; y2 = ig.bottom_thk + ig.contact_thk;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_busbar_base', ui);
    ig.cu_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];
end

% --- Baseplate rectangle ---
ig.xmin = min(ig.x_centers) - (ig.a + ig.busbar_thk) - ig.margin;
ig.xmax = max(ig.x_centers) + (ig.a + ig.busbar_thk) + ig.margin;
R = [3;4; ig.xmin; ig.xmax; ig.xmax; ig.xmin; -ig.base_thk; -ig.base_thk; 0; 0];
gd(:,end+1) = R; idx = idx+1; 
names{end+1} = 'baseplate';
ig.al_centroid(end+1,:) = [(ig.xmin+ig.xmax)/2, (-ig.base_thk)/2];

% --- Heater/contact strips ---
y1 = 0; y2 = ig.contact_thk;
for ui = 1:ig.numU
    xc = ig.x_centers(ui);
    x1 = xc - ig.a; x2 = xc + ig.a;
    R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2]; 
    gd(:,end+1) = R; idx = idx+1; 
    names{end+1} = sprintf('U%d_thermal', ui);
end

% --- Cells inside U channels ---
cell_w = ig.cell_length;  % width in x-direction
cell_h = ig.cell_diameter;  % height in y-direction

for ui = 1:ig.numU
    xc = ig.x_centers(ui);  % center x of this U
    yc_centers = (0.5:1:ig.Nrows) .* ((ig.module_h - ig.base_thk) / (ig.Nrows)) + ig.base_thk + ig.contact_thk;  % same vertical spacing as heaters

    for r = 1:ig.Nrows
        yc = yc_centers(r);   % center y position of this row
        
        % rectangle spanning from xc-a → xc+a (inner U walls)
        x1 = xc - cell_w/2; 
        x2 = xc + cell_w/2;
        y1 = yc - cell_h/2;
        y2 = yc + cell_h/2;
        
        R = [3;4; x1; x2; x2; x1; y1; y1; y2; y2];
        gd(:,end+1) = R;
        idx = idx+1;
        names{end+1} = sprintf('U%d_R%d_cell',ui,r);
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
        ig.al_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];
    end
end



% --- Store centroids ---
% --- TIM/heater strips ---
for ui = 1:ig.numU
    xc = ig.x_centers(ui);
    y1 = 0; y2 = ig.contact_thk;
    xc = ig.x_centers(ui);
    x1 = xc - ig.a; x2 = xc + ig.a;
    ig.tim_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];
end

% --- Cells ---
cell_w = ig.cell_length;
cell_h = ig.cell_diameter;
for ui = 1:ig.numU
    xc = ig.x_centers(ui);
    yc_centers = (0.5:1:ig.Nrows) .* ((ig.module_h - ig.base_thk) / ig.Nrows) + ig.base_thk;
    for r = 1:ig.Nrows
        yc = yc_centers(r);
        x1 = xc - cell_w/2; x2 = xc + cell_w/2;
        y1 = yc - cell_h/2; y2 = yc + cell_h/2;
        ig.cell_centroid(end+1,:) = [(x1+x2)/2, (y1+y2)/2];
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
