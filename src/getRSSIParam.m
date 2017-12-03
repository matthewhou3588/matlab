function [Z, K] = getRSSIParam(x, y)
%% use least square method in offline to get Z and K in 
%% equation RSSI = Z*log(w)+K for every access point 
%% (called beacon node in paper)

% x = [1 2 3 4 5 6 7];
% y = [5 7 9 11 13 15 17];
% x = 0:0.1:1;
% y = 2*x+3;
[P, S, mu] = polyfit(x, y, 1);
%z = polyval(P, x)
Z = P(1,1);
K = P(1,2);


