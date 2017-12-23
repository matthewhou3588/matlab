function [ unknown_node_location ] = online_stage( env, grid, unknown_node )
% this online_stage is for locate one unknown node


%% 1.
%% RSSI compute in online
%% the first step: beacon node is grid
%  the RSSI values received by the beacon are computed. Only the lowest values
%  (in absolute value), under the 25th percentile are considered.
rowpercentile = quantile(grid.rssi, 0.25, 2);  % for every beacon node, get the 25th percentile of RSSI


%% 2. 
%% second step, the distance between the unknown node and the beacon n is estimated
anchor_unknown_estimated_dist = zeros(1, grid.num);
for k = 1:grid.num
    anchor_unknown_estimated_dist(1,k) = 10^3*sqrt(sum((repmat(grid.points(k,:),1,1)-unknown_node.points).^2,2)).';  % the distance from the kth grid to all access points
end

%% 3.
%% third step 
% for each cell, the value w藛n 鈭抴n(i, j) is
% associated to it in order to create an error map related to anchor n
s = 0.1;
startloc = -sqrt(env.area)/2 + s/2;
endloc = sqrt(env.area)/2 - s/2;

xpos = linspace(startloc, endloc, sqrt(env.area)/s);
ypos = linspace(startloc, endloc, sqrt(env.area)/s);

cell_center = zeros(sqrt(env.area)/s * sqrt(env.area)/s, 2);  % the centers of map cells
for i=1:length(xpos)
    for j=1:length(ypos)
        cell_center((i-1)*length(xpos)+j, 1) = xpos(1,i);
        cell_center((i-1)*length(xpos)+j, 2) = ypos(1,j);
    end
end 
% compute wn(i, j)
cellcenter_anchor_dist = zeors(sqrt(env.area)/s * sqrt(env.area)/s, grid.num); % every row : the distance of center of cell i to every grid node(anrchor node)
for k = 1:size(cellcenter_anchor_dist, 1)  
    cellcenter_anchor_dist(k,:) = 10^3*sqrt(sum((repmat(cell_center(k,:), grid.num, 1) - grid.points).^2, 2)).'; 
end

err_map_dist = zeors(sqrt(env.area)/s * sqrt(env.area)/s, grid.num);
for k = 1:size(err_map_dist, 1)  
    % err_map_dist(k,:) = abs(anchor_unknown_estimated_dist - cellcenter_anchor_dist(k,:));
    err_map_dist(k,:) = (anchor_unknown_estimated_dist - cellcenter_anchor_dist(k,:)).^2;
end


%% 4.
%% weighting every anchor node's reliability, 
%% the weighting result is the outpout of the flc2

% 暂时假设所有的grid的reliability如下：
anchor_reliability = unifrnd(0, 1, 1, grid.num); 

%% 5.
W_i_j = zeros(1, size(err_map_dist, 1));  % every cell's result, and then find the min from all the cell
for k = 1:size(err_map_dist, 1)  
    tmp = 0;
    for j = 1:grid.num;      %% tmp is the anchor node j's W
        tmp = tmp + anchor_reliability(1,j) * err_map_dist(k,j);  
    end
    W_i_j(1, k) = tmp;  % the kth cell's W
end

% find the min W (the least error) from all cell
[min_error, cell_index] = min(W_i_j);
unknown_node_location = cell_center(cell_index,:);


end
