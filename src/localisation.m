clear
close all 

%% variables that you could change
ap.num = 36;    % the number of access points (must be square of interger)
grid.num = 900;  % the number of fingerprinting grids (must be square of interger)
ue.density = ap.num;  % the density of user equipments
ap.tx_power = 0.2*ones(1,ap.num);  % the tx power of access points in Watts

%% constants that you do not need to change
env.area = 1;   % the playground area in square km
env.dc = [6,4*6*1.25/(3*10^8/(2*10^9))];  % multi-slope pathloss threshold

%% generating the map
ap.xpos = linspace(-sqrt(env.area)/2,sqrt(env.area)/2,2*sqrt(ap.num)-1+2);
ap.ypos = linspace(-sqrt(env.area)/2,sqrt(env.area)/2,2*sqrt(ap.num)-1+2);
ap.xpos = ap.xpos(2:2:end);
ap.ypos = ap.ypos(2:2:end);
ap.points = [repmat(ap.xpos,1,sqrt(ap.num));reshape(repmat(ap.ypos,sqrt(ap.num),1),[ap.num,1]).'].';

grid.xpos = linspace(-sqrt(env.area)/2,sqrt(env.area)/2,2*sqrt(grid.num)-1+2);
grid.ypos = linspace(-sqrt(env.area)/2,sqrt(env.area)/2,2*sqrt(grid.num)-1+2);
grid.xpos = grid.xpos(2:2:end);
grid.ypos = grid.ypos(2:2:end);
grid.points = [repmat(grid.xpos,1,sqrt(grid.num));reshape(repmat(grid.ypos,sqrt(grid.num),1),[grid.num,1]).'].';

ue.num = poissrnd(ue.density*env.area);  % the number of user equipments
ue.points = unifrnd(-sqrt(env.area)/2,sqrt(env.area)/2,ue.num,2);   % the position of user equipments

% figure
% plot(ap.points(:,1),ap.points(:,2),'b^','MarkerSize',4)
% hold on
% plot(grid.points(:,1),grid.points(:,2),'ks','MarkerSize',4)
% axis([-sqrt(env.area)/2,sqrt(env.area)/2,-sqrt(env.area)/2,sqrt(env.area)/2],'square')
% legend('AP','grid');

% figure
% plot(ap.points(:,1),ap.points(:,2),'b^','MarkerSize',4)
% hold on
% plot(ue.points(1,1),ue.points(1,2),'ro','MarkerSize',4)
% axis([-sqrt(env.area)/2,sqrt(env.area)/2,-sqrt(env.area)/2,sqrt(env.area)/2],'square')
% legend('AP','UE'); 


%% generating distance, pathloss and rssi
for k = 1:grid.num
    grid.dist(k,:) = 10^3*sqrt(sum((repmat(grid.points(k,:),ap.num,1)-ap.points).^2,2)).';  % the distance from the kth grid to all access points
    grid.pathloss(k,grid.dist(k,:) <= env.dc(1)) = grid.dist(k,grid.dist(k,:) <= env.dc(1)).^(-0);  % the pathloss from the kth grid to all access points
    grid.pathloss(k,grid.dist(k,:) > env.dc(1) & grid.dist(k,:) <= env.dc(2)) = env.dc(1)^2*grid.dist(k,grid.dist(k,:) > env.dc(1) & grid.dist(k,:) <= env.dc(2)).^(-2);
    grid.pathloss(k,grid.dist(k,:) > env.dc(2)) = env.dc(1)^2*env.dc(2)^2*grid.dist(k,grid.dist(k,:) > env.dc(2)).^(-4);
    grid.rssi(k,:) = grid.pathloss(k,:).*ap.tx_power;    % the rssi from the kth grid to all access points
end

for k = 1:ue.num
    ue.dist(k,:) = 10^3*sqrt(sum((repmat(ue.points(k,:),ap.num,1)-ap.points).^2,2)).';  % the distance from the kth user equipment to all access points
    ue.pathloss(k,ue.dist(k,:) <= env.dc(1)) = ue.dist(k,ue.dist(k,:) <= env.dc(1)).^(-0);  % the pathloss from the kth user equipment to all access points
    ue.pathloss(k,ue.dist(k,:) > env.dc(1) & ue.dist(k,:) <= env.dc(2)) = env.dc(1)^2*ue.dist(k,ue.dist(k,:) > env.dc(1) & ue.dist(k,:) <= env.dc(2)).^(-2);
    ue.pathloss(k,ue.dist(k,:) > env.dc(2)) = env.dc(1)^2*env.dc(2)^2*ue.dist(k,ue.dist(k,:) > env.dc(2)).^(-4);
    ue.rssi(k,:) = ue.pathloss(k,:).*ap.tx_power;    % the rssi from the kth user equipment to all access points    
end

%% clear 
% clear env;
% clear ;


    for n=1:ue.num 
       for i=1:grid.num
           error(i,:) = abs(ue.rssi(n,:)-grid.rssi(i,:));
           error_sum(i,:) = sum(error(i,:)); %Signal Processing
       end
    [min_error_sum(n,:),temp_position_index] = min(error_sum); % closest point    
    position_index(1,n)=temp_position_index;
    end
    
    for j=1:ue.num
    
            location(j,1) = grid.points(position_index(1,j),1);
            location(j,2) = grid.points(position_index(1,j),2);
    end

figure
plot(location(:,1),location(:,2),'bo','MarkerSize',4)
hold on
plot(ue.points(:,1),ue.points(:,2),'ro','MarkerSize',4)
axis([-sqrt(env.area)/2,sqrt(env.area)/2,-sqrt(env.area)/2,sqrt(env.area)/2],'square')
legend('localization','ue')
 

