% track_generator.m
% Generates a track from a centreline for LTS model
% Original: Lewis Blake
% Edits: use all SVG points + dedupe near-duplicates to avoid NaNs/0-length
%        and minor plotting handle fix.

%function [dels,curv,track_length] = track_generator()

close all
clear

halfwidth = 1.5;                 % track halfwidth, metres
carhalfwidth = halfwidth/1.5;
L_max = carhalfwidth;
L_min = -carhalfwidth;

mat = loadsvg("FSUK2025.svg", 1, false);  % TRACK INPUT

% --- Build centreline from *all* SVG points with a tiny dedupe guard -----
xy_raw = mat{1,1};                        % Nx2
if size(xy_raw,2) ~= 2
    error('Expected mat{1,1} to be Nx2 coordinates.');
end

% Deduplicate consecutive points that are effectively identical
d = vecnorm(diff(xy_raw), 2, 2);
nz = d(d > 0);
if isempty(nz)
    medseg = 0;
else
    medseg = median(nz);
end
eps_close = max(1e-9, 1e-6*medseg);       % scale-aware tolerance

keep = [true; d > eps_close];
xy = xy_raw(keep,:);

% If path is closed and last repeats first, drop the duplicate closing point
if norm(xy(end,:) - xy(1,:)) < eps_close
    xy(end,:) = [];
end

if size(xy,1) < 3
    error('Not enough unique points after dedupe.');
end

% Allocate track (keep 3 cols initially so length(track) is #rows)
track = zeros(size(xy,1), 3);
track(:,1:2) = xy;
track_length = 0;

% -------------------- Geom precompute: normals, limits --------------------
for i = 2:(length(track)-1)
    % arc step and running length
    track(i,8) = sqrt((track(i+1,2)-track(i,2))^2 + (track(i+1,1)-track(i,1))^2);
    track_length = track_length + track(i,8);

    % local perpendicular (from i-1 to i+1)
    A = [track(i+1,1), track(i+1,2)];
    B = [track(i-1,1), track(i-1,2)];
    AB = B - A;
    AB = AB / norm(AB);
    ABperp = AB*[0 -1; 1 0];

    track(i,13) = ABperp(1);
    track(i,14) = ABperp(2);

    % track limits (full width)
    track(i,4) = track(i,1) + halfwidth*ABperp(1);
    track(i,5) = track(i,2) + halfwidth*ABperp(2);
    track(i,6) = track(i,1) - halfwidth*ABperp(1);
    track(i,7) = track(i,2) - halfwidth*ABperp(2);

    % car envelope limits (narrower)
    track(i,9)  = track(i,1) + carhalfwidth*ABperp(1);
    track(i,10) = track(i,2) + carhalfwidth*ABperp(2);
    track(i,11) = track(i,1) - carhalfwidth*ABperp(1);
    track(i,12) = track(i,2) - carhalfwidth*ABperp(2);
end

L = zeros(length(track),1);
L_temp = zeros(length(track),1);
race_line = [track(:,1), track(:,2)];

% removes null island points (first/last sentinel row)
track2 = track(2:(end-1),:);
race_line2 = race_line(2:(end-1),:); %#ok<NASGU>

% ------------------------------ Plotting ----------------------------------
figure; hold on; axis equal
h_center = plot(track2(:,1), track2(:,2), 'r-');  % centreline
set(h_center, 'LineWidth', 2);

plot(track2(:,4), track2(:,5), 'k-', 'LineWidth', 2); % left limit
plot(track2(:,6), track2(:,7), 'k-', 'LineWidth', 2); % right limit
title('Track and Limits');
xlabel('x [m]'); ylabel('y [m]');

% --------------------------- Optimisation ---------------------------------
curv_param = 4;
dist_param = 2;
threshold = 5e-3;
count = 0;
check_count = 0;
testcount = 0; %#ok<NASGU>
delL = 10;

while delL > threshold
    L_old = L;

    for j = 2:(length(track)-1)
        for i = [17,33,55,75,155]    % gaussian kernel widths (odd)
            % pre-alloc (length(track) = #points)
            theta2_temp = zeros(length(track),1);
            theta1_temp = zeros(length(track),1);
            deltheta_temp = zeros(length(track),1);
            s2_temp = zeros(length(track),1);
            s1_temp = zeros(length(track),1);
            dels_temp = zeros(length(track),1);
            curv_temp = zeros(length(track),1);
            theta2 = zeros(length(track),1);
            theta1 = zeros(length(track),1);
            deltheta = zeros(length(track),1);
            s2 = zeros(length(track),1);
            s1 = zeros(length(track),1);
            dels = zeros(length(track),1);
            curv = zeros(length(track),1);
            costsum_temp = 0;
            costsum = 0;

            track_temp = [race_line(:,1), race_line(:,2)];
            L_temp = L;

            [x, y] = gauss_gen(i);   % indices x, weights y (length i), i must be odd

            % positive shift along local normal
            for num = x
                if ((num+j) > 0) && ((num+j) < length(track))
                    track_temp(num+j,1) = track_temp(num+j,1) + y(num+((i+1)/2))*track(num+j,13);
                    track_temp(num+j,2) = track_temp(num+j,2) + y(num+((i+1)/2))*track(num+j,14);
                    L_temp(num+j)       = L_temp(num+j)       + y(num+((i+1)/2));
                end
            end

            % compute costs (temp vs current)
            for p = 2:(length(track)-1)
                % curvature of temp track
                theta2_temp(p,1) = atand((track_temp(p+1,2)-track_temp(p,2))/((track_temp(p+1,1)-track_temp(p,1))));
                theta1_temp(p,1) = atand((track_temp(p,2)-track_temp(p-1,2))/((track_temp(p,1)-track_temp(p-1,1))));
                deltheta_temp(p,1) = theta2_temp(p,1) - theta1_temp(p,1);
                if deltheta_temp(p,1) > 90,  deltheta_temp(p,1) = deltheta_temp(p,1) - 180; end
                if deltheta_temp(p,1) < -90, deltheta_temp(p,1) = deltheta_temp(p,1) + 180; end
                s2_temp(p,1) = sqrt((track_temp(p+1,1)-track_temp(p,1))^2 + (track_temp(p+1,2)-track_temp(p,2))^2);
                s1_temp(p,1) = sqrt((track_temp(p-1,1)-track_temp(p,1))^2 + (track_temp(p-1,2)-track_temp(p,2))^2);
                dels_temp(p,1) = (s1_temp(p,1)+s2_temp(p,1))/2;
                curv_temp(p,1) = (deltheta_temp(p,1)) / (dels_temp(p,1));
                costsum_temp    = costsum_temp + (curv_temp(p,1)^curv_param) * (dels_temp(p,1)^dist_param);

                % curvature of current race_line
                theta2(p,1) = atand((race_line(p+1,2)-race_line(p,2))/((race_line(p+1,1)-race_line(p,1))));
                theta1(p,1) = atand((race_line(p,2)-race_line(p-1,2))/((race_line(p,1)-race_line(p-1,1))));
                deltheta(p,1) = theta2(p,1) - theta1(p,1);
                if deltheta(p,1) > 90,  deltheta(p,1) = deltheta(p,1) - 180; end
                if deltheta(p,1) < -90, deltheta(p,1) = deltheta(p,1) + 180; end
                s2(p,1) = sqrt((race_line(p+1,1)-race_line(p,1))^2 + (race_line(p+1,2)-race_line(p,2))^2);
                s1(p,1) = sqrt((race_line(p-1,1)-race_line(p,1))^2 + (race_line(p-1,2)-race_line(p,2))^2);
                dels(p,1) = (s1(p,1)+s2(p,1))/2;
                curv(p,1) = (deltheta(p,1)) / (dels(p,1));
                costsum    = costsum    + (curv(p,1)^curv_param) * (dels(p,1)^dist_param);
            end

            check_count = check_count + 1;
            if (costsum_temp < costsum) && (max(L_temp) < L_max)
                race_line = track_temp;
                L = L_temp;
                count = count + 1;
            end

            % negative shift along local normal
            theta2_temp(:) = 0; theta1_temp(:) = 0; deltheta_temp(:) = 0;
            s2_temp(:) = 0; s1_temp(:) = 0; dels_temp(:) = 0; curv_temp(:) = 0;
            theta2(:) = 0; theta1(:) = 0; deltheta(:) = 0; s2(:) = 0; s1(:) = 0; dels(:) = 0; curv(:) = 0;
            costsum_temp = 0; costsum = 0;
            track_temp = [race_line(:,1), race_line(:,2)];
            L_temp = L;

            for num = x
                if (num+j > 0) && (num+j < length(track))
                    track_temp(num+j,1) = track_temp(num+j,1) - y(num+((i+1)/2))*track(num+j,13);
                    track_temp(num+j,2) = track_temp(num+j,2) - y(num+((i+1)/2))*track(num+j,14);
                    L_temp(num+j)       = L_temp(num+j)       - y(num+((i+1)/2));
                end
            end

            for p = 2:(length(track)-1)
                % curvature of temp track
                theta2_temp(p,1) = atand((track_temp(p+1,2)-track_temp(p,2))/((track_temp(p+1,1)-track_temp(p,1))));
                theta1_temp(p,1) = atand((track_temp(p,2)-track_temp(p-1,2))/((track_temp(p,1)-track_temp(p-1,1))));
                deltheta_temp(p,1) = theta2_temp(p,1) - theta1_temp(p,1);
                if deltheta_temp(p,1) > 90,  deltheta_temp(p,1) = deltheta_temp(p,1) - 180; end
                if deltheta_temp(p,1) < -90, deltheta_temp(p,1) = deltheta_temp(p,1) + 180; end
                s2_temp(p,1) = sqrt((track_temp(p+1,1)-track_temp(p,1))^2 + (track_temp(p+1,2)-track_temp(p,2))^2);
                s1_temp(p,1) = sqrt((track_temp(p-1,1)-track_temp(p,1))^2 + (track_temp(p-1,2)-track_temp(p,2))^2);
                dels_temp(p,1) = (s1_temp(p,1)+s2_temp(p,1))/2;
                curv_temp(p,1) = (deltheta_temp(p,1)) / (dels_temp(p,1));
                costsum_temp    = costsum_temp + (curv_temp(p,1)^curv_param) * (dels_temp(p,1)^dist_param);

                % curvature of current race_line
                theta2(p,1) = atand((race_line(p+1,2)-race_line(p,2))/((race_line(p+1,1)-race_line(p,1))));
                theta1(p,1) = atand((race_line(p,2)-race_line(p-1,2))/((race_line(p,1)-race_line(p-1,1))));
                deltheta(p,1) = theta2(p,1) - theta1(p,1);
                if deltheta(p,1) > 90,  deltheta(p,1) = deltheta(p,1) - 180; end
                if deltheta(p,1) < -90, deltheta(p,1) = deltheta(p,1) + 180; end
                s2(p,1) = sqrt((race_line(p+1,1)-race_line(p,1))^2 + (race_line(p+1,2)-race_line(p,2))^2);
                s1(p,1) = sqrt((race_line(p-1,1)-race_line(p,1))^2 + (race_line(p-1,2)-race_line(p,2))^2);
                dels(p,1) = (s1(p,1)+s2(p,1))/2;
                curv(p,1) = (deltheta(p,1)) / (dels(p,1));
                costsum    = costsum    + (curv(p,1)^curv_param) * (dels(p,1)^dist_param);
            end

            check_count = check_count + 1;
            if (costsum_temp < costsum) && (min(L_temp) > L_min)
                race_line = track_temp;
                L = L_temp;
                count = count + 1;
            end
        end
    end

    delL = abs(norm(L) - norm(L_old));
    % fprintf('Outer iter: delL=%.6g, accepted=%d\n', delL, count); % optional
end

% --------------------------- Final plotting --------------------------------
plot(race_line(:,1), race_line(:,2), 'b-', 'LineWidth', 1.5)
plot(track2(:,9),  track2(:,10), 'k--')
plot(track2(:,11), track2(:,12), 'k--')

print -depsc endurancetrack

% --------------------------- Output & save ---------------------------------
curv = smooth(curv, 10);
save('curve.mat','curv', 'dels','track_length');

%end
