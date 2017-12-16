clear 
close all 

%% variables that you could change
ap.num = 25;    % the number of access points (must be square of interger)
grid.num = 100;  % the number of fingerprinting grids (must be square of interger)
ue.density = ap.num;  % the density of user equipments
ap.tx_power = 0.2*ones(1,ap.num);  % the tx power of access points in Watts
% receiver.num = 1; % the number of the recevier node
% receiver.pos = [0 0]; % the location of the recevier node



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
% 
% figure
% plot(ap.points(:,1),ap.points(:,2),'b^','MarkerSize',4)
% hold on
% plot(ue.points(:,1),ue.points(:,2),'ro','MarkerSize',4)
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

%% calculte the distance, rssi between every access point and the receiver node
% for k = 1:ap.num
%     receiver.dist(k,:) = 10^3*sqrt(sum(ap.points(k,:).^2)).';  % the distance from the kth access point to receiver node
%     ue.pathloss(k,ue.dist(k,:) <= env.dc(1)) = ue.dist(k,ue.dist(k,:) <= env.dc(1)).^(-0);  % the pathloss from the kth access point to receiver node
%     ue.pathloss(k,ue.dist(k,:) > env.dc(1) & ue.dist(k,:) <= env.dc(2)) = env.dc(1)^2*ue.dist(k,ue.dist(k,:) > env.dc(1) & ue.dist(k,:) <= env.dc(2)).^(-2);
%     ue.pathloss(k,ue.dist(k,:) > env.dc(2)) = env.dc(1)^2*env.dc(2)^2*ue.dist(k,ue.dist(k,:) > env.dc(2)).^(-4);
%     ue.rssi(k,:) = ue.pathloss(k,:).*ap.tx_power;    % the rssi from the kth access point to receiver node
% end


%% clear 
% clear env;
% clear k;

%% get z and k for every ap
Z = zeros(1, ap.num);
K = zeros(1, ap.num);
for k = 1:ap.num
    x = log10(grid.dist(:,k).');  
    y = grid.rssi(:,k).';
    [Z(1,k), K(1,k)] = getRSSIParam(x, y);
end

Z
K

%% get every ap's reliability in offline
% [ result, z_output, k_output ] = flc1( Z, K );


%% get every ap's membership function of flc1 in offline
membership_args = ones(ap.num, 9);
for k = 1:2  %ap.num
    retval = run_pso(Z(1,k), K(1,k));
    
    membership_args(k,1) = 0;                % aL
    membership_args(k,2) = retval(1,2);  % bL
    membership_args(k,3) = retval(1,1);  % cL

    membership_args(k,4) = retval(1,3);  % aM
    membership_args(k,5) = 50;               % bM
    membership_args(k,6) = retval(1,4);  % cM

    membership_args(k,7) = retval(1,5);  % aH
    membership_args(k,8) = retval(1,6);  % bH
    membership_args(k,9) = 100;              % cH
    
    % print
    membership_args(k,:)
end

%% RSSI compute in online
%% the first step: beacon node is grid
%  the RSSI values received by the beacon are computed. Only the lowest values
%  (in absolute value), under the 25th percentile are considered.
rowpercentile = quantile(grid.rssi, 0.25, 2);  % for every beacon node, get the 25th percentile of all ap

%% second step, the distance between the unknown node and the beacon n is estimated
w_unknow_beancon = ue.dist;

%% third step 
% for each cell, the value kwˆn −wn(i, j)k is
% associated to it in order to create an error map related to anchor n

for i=1:ue.num  % loccate every unknown node each for. here we think unknow node(paper) is ue
    
   errmap = zeros(grid.num, ap.num);   % the ith ue's error map between this ue and every grid cell
   for j=1:grid.num
       errmap(j,:) = w_unknow_beancon(i,:) - grid.dist(j,:); 
   end
    
   %% final step
   proximity_index_matrix = zeros(grid.num, ap.num); 
   strongestRSSIfromUnknow = min(abs(ue.rssi(i,:)));
   for j=1:grid.num      
       for t=1:ap.num
           proximity_index_matrix(j,t) = strongestRSSIfromUnknow / grid.rssi(j,t);
       end       
   end

   % normalize proximity_index
   maxv = max(proximity_index_matrix);
   minv = min(proximity_index_matrix);   
   proximity_index_matrix_normalized = (proximity_index_matrix - minv) / (maxv - minv);     

   
   
end


