function [Z, K] = getRSSIParam(x, y)
%% use least square method in offline to get Z and K in 
%% equation RSSI = Z*log(w)+K for every 
%% (called beacon node in paper)

%原始数据  
% x=[163     123     150      123     141];  
% y=[186     126     172      125     148];  
% [P, S, mu] = polyfit(x, y, 1);
% P
% Z = P(1,1);
% K = P(1,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x);
x2=sum(x.^2);       % 求Σ(xi^2)  
x1=sum(x);          % 求Σ(xi)  
x1y1=sum(x.*y);     % 求Σ(xi*yi)  
y1=sum(y);          % 求Σ(yi)  
  
a=(n*x1y1-x1*y1)/(n*x2-x1*x1);      %解出直线斜率b=(y1-a*x1)/n  
b=(y1-a*x1)/n;                      %解出直线截距  

Z = a;
K = b;

% %作图  
% % 先把原始数据点用蓝色十字描出来  
% figure  
% plot(x,y,'+');        
% hold on  
% % 用红色绘制拟合出的直线  
% px=linspace(min(x),max(x));%这里直线区间根据自己实际需求改写  
% py=a*px+b;  
% plot(px,py,'r');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %作图  
% % 先把原始数据点用蓝色十字描出来  
% figure  
% plot(x,y,'+');        
% hold on  
% % 用红色绘制拟合出的直线  
% px=linspace(min(x),max(x),45);%这里直线区间根据自己实际需求改写  
% py=Z*px+K;  
% plot(px,py,'r');  


