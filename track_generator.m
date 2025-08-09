%track_generator.m
%Generates a track from a centreline for LTS model
%Lewis Blake

%function [dels,curv,track_length] = track_generator()

close all
clear

halfwidth = 1.5; %track halfwidth, metres
carhalfwidth = halfwidth/1.5;
L_max = carhalfwidth;
L_min = -carhalfwidth;
mat = loadsvg("FSUK2025.svg",1,false); %TRACK INPUT
track = zeros(length(mat{1,1})/2,3);
track_length = 0;

for i = 1:2:length(mat{1,1})  %USE THIS FOR ENDURANCE TRACK
    track(((i+1)/2),1) = mat{1,1}(i,1);
    track(((i+1)/2),2) = mat{1,1}(i,2);
end




for i = 2:(length(track)-1)

    track(i,8) = sqrt((track(i+1,2)-track(i,2))^2+(track(i+1,1)-track(i,1))^2);
    track_length = track_length + track(i,8);

    A = [track(i+1,1),track(i+1,2)];
    B = [track(i-1,1),track(i-1,2)];
    AB = B - A;
    AB = AB/norm(AB);
    ABperp = AB*[0 -1;1 0];
    track(i,13) = ABperp(1);
    track(i,14) = ABperp(2);

    track(i,4) = track(i,1) + halfwidth*ABperp(1);
    track(i,5) = track(i,2) + halfwidth*ABperp(2);
    track(i,6) = track(i,1) - halfwidth*ABperp(1);
    track(i,7) = track(i,2) - halfwidth*ABperp(2);

    track(i,9) = track(i,1) + carhalfwidth*ABperp(1);
    track(i,10) = track(i,2) + carhalfwidth*ABperp(2);
    track(i,11) = track(i,1) - carhalfwidth*ABperp(1);
    track(i,12) = track(i,2) - carhalfwidth*ABperp(2);



end


L = zeros(length(track),1);
L_temp = zeros(length(track),1);
race_line = [track(:,1),track(:,2)];



%removes null island points
track2 = track(2:(end-1),:);
race_line2 = race_line(2:(end-1),:);


%plots centreline
plot(track2(:,1),track2(:,2),"red-")
hold on
%plot(race_line(:,1),race_line(:,2),"blue-")
%plots track limits
LineHandle.NodeChildren.LineWidth = 2;
plot(track2(:,4),track2(:,5),"black-")
plot(track2(:,6),track2(:,7),"black-")

% marks widest line for racing line
% plot(track2(:,9),track2(:,10),"black--")
% plot(track2(:,11),track2(:,12),"black--")


%PLOTTING MORE LINES

curv_param = 4;
dist_param = 2;
threshold = 5*10^-3;
count = 0;
check_count = 0;
testcount = 0;
delL = 10;
while delL > threshold
    L_old = L;
    %for i = [5,15,25,55]
    for j = 2:(length(track)-1)
        for i = [17,33,55,75,155]
        %for j = 2:(length(track)-1)

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
            track_temp = [race_line(:,1),race_line(:,2)];
            L_temp = L;
            [x,y] = gauss_gen(i);
            for num = x
                if ((num+j) > 0) && (num+j < length(track))
                    track_temp(num+j,1) = track_temp(num+j,1)+y(num+((i+1)/2))*track(num+j,13);
                    track_temp(num+j,2) = track_temp(num+j,2)+y(num+((i+1)/2))*track(num+j,14);
                    L_temp(num+j) = L_temp(num+j)+y(num+((i+1)/2));
                    
                end
            end
            for p = 2:(length(track)-1)
                %finding curvature of temp track
                theta2_temp(p,1) = atand((track_temp(p+1,2)-track_temp(p,2))/((track_temp(p+1,1)-track_temp(p,1))));
                theta1_temp(p,1) = atand((track_temp(p,2)-track_temp(p-1,2))/((track_temp(p,1)-track_temp(p-1,1))));
                deltheta_temp(p,1) = theta2_temp(p,1)-theta1_temp(p,1);
                if deltheta_temp(p,1) > 90
                    deltheta_temp(p,1)=deltheta_temp(p,1)-180;
                end
                if deltheta_temp(p,1) < -90
                    deltheta_temp(p,1)=deltheta_temp(p,1)+180;
                end
                s2_temp(p,1) = sqrt((track_temp(p+1,1)-track_temp(p,1))^2+(track_temp(p+1,2)-track_temp(p,2))^2);
                s1_temp(p,1) = sqrt((track_temp(p-1,1)-track_temp(p,1))^2+(track_temp(p-1,2)-track_temp(p,2))^2);
                dels_temp(p,1) = (s1_temp(p,1)+s2_temp(p,1))/2;
                curv_temp(p,1) = (deltheta_temp(p,1))/(dels_temp(p,1));
                costsum_temp = costsum_temp + (curv_temp(p,1)^curv_param)*(dels_temp(p,1)^dist_param);
                
                
                %finding curvature
                theta2(p,1) = atand((race_line(p+1,2)-race_line(p,2))/((race_line(p+1,1)-race_line(p,1))));
                theta1(p,1) = atand((race_line(p,2)-race_line(p-1,2))/((race_line(p,1)-race_line(p-1,1))));
                deltheta(p,1) = theta2(p,1)-theta1(p,1);
                if deltheta(p,1) > 90
                    deltheta(p,1)=deltheta(p,1)-180;
                end
                if deltheta(p,1) < -90
                    deltheta(p,1)=deltheta(p,1)+180;
                end
                s2(p,1) = sqrt((race_line(p+1,1)-race_line(p,1))^2+(race_line(p+1,2)-race_line(p,2))^2);
                s1(p,1) = sqrt((race_line(p-1,1)-race_line(p,1))^2+(race_line(p-1,2)-race_line(p,2))^2);
                dels(p,1) = (s1(p,1)+s2(p,1))/2;
                curv(p,1) = (deltheta(p,1))/(dels(p,1));
                costsum = costsum + (curv(p,1)^curv_param)*(dels(p,1)^dist_param);
            end
            check_count = check_count+1;
            if (costsum_temp < costsum) && (max(L_temp) < L_max)
                race_line = track_temp;
                L = L_temp;
                
                count=count+1;
                
            end
            %L_old = L;
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
            track_temp = [race_line(:,1),race_line(:,2)];
            L_temp = L;
            for num = x
                if (num+j > 0) && (num+j < length(track))
                    track_temp(num+j,1) = track_temp(num+j,1)-y(num+((i+1)/2))*track(num+j,13);
                    track_temp(num+j,2) = track_temp(num+j,2)-y(num+((i+1)/2))*track(num+j,14);
                    L_temp(num+j) = L_temp(num+j) - y(num+((i+1)/2));
                end
            end
            for p = 2:(length(track)-1)
                %finding curvature of temp track
                theta2_temp(p,1) = atand((track_temp(p+1,2)-track_temp(p,2))/((track_temp(p+1,1)-track_temp(p,1))));
                theta1_temp(p,1) = atand((track_temp(p,2)-track_temp(p-1,2))/((track_temp(p,1)-track_temp(p-1,1))));
                deltheta_temp(p,1) = theta2_temp(p,1)-theta1_temp(p,1);
                if deltheta_temp(p,1) > 90
                    deltheta_temp(p,1)=deltheta_temp(p,1)-180;
                end
                if deltheta_temp(p,1) < -90
                    deltheta_temp(p,1)=deltheta_temp(p,1)+180;
                end
                s2_temp(p,1) = sqrt((track_temp(p+1,1)-track_temp(p,1))^2+(track_temp(p+1,2)-track_temp(p,2))^2);
                s1_temp(p,1) = sqrt((track_temp(p-1,1)-track_temp(p,1))^2+(track_temp(p-1,2)-track_temp(p,2))^2);
                dels_temp(p,1) = (s1_temp(p,1)+s2_temp(p,1))/2;
                curv_temp(p,1) = (deltheta_temp(p,1))/(dels_temp(p,1));
                costsum_temp = costsum_temp + (curv_temp(p,1)^curv_param)*(dels_temp(p,1)^dist_param);
                
                
                %finding curvature
                theta2(p,1) = atand((race_line(p+1,2)-race_line(p,2))/((race_line(p+1,1)-race_line(p,1))));
                theta1(p,1) = atand((race_line(p,2)-race_line(p-1,2))/((race_line(p,1)-race_line(p-1,1))));
                deltheta(p,1) = theta2(p,1)-theta1(p,1);
                if deltheta(p,1) > 90
                    deltheta(p,1)=deltheta(p,1)-180;
                end
                if deltheta(p,1) < -90
                    deltheta(p,1)=deltheta(p,1)+180;
                end
                s2(p,1) = sqrt((race_line(p+1,1)-race_line(p,1))^2+(race_line(p+1,2)-race_line(p,2))^2);
                s1(p,1) = sqrt((race_line(p-1,1)-race_line(p,1))^2+(race_line(p-1,2)-race_line(p,2))^2);
                dels(p,1) = (s1(p,1)+s2(p,1))/2;
                curv(p,1) = (deltheta(p,1))/(dels(p,1));
                costsum = costsum + (curv(p,1)^curv_param)*(dels(p,1)^dist_param);
                
            end
            check_count = check_count+1;
            if (costsum_temp < costsum) && (min(L_temp) > L_min)
                race_line = track_temp;
                L = L_temp;
                count=count+1;
                
            end            
        end   
    end
        
    delL = abs(norm(L)-norm(L_old))

end



axis equal
plot(race_line(:,1),race_line(:,2),"blue-")

%marks widest line for racing line
plot(track2(:,9),track2(:,10),"black--")
plot(track2(:,11),track2(:,12),"black--")


print -depsc endurancetrack

curv = smooth(curv,10);

save('curve.mat','curv', 'dels','track_length');

%end