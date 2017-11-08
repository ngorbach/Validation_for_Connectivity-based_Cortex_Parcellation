function [x,t] = OU_process(starting_point,target_point,final_time,sig)

th = 1;
dt = 1e-2;
t = 0:dt:final_time;             % Time vector
%rng(1);                 % Set random seed
W = zeros(1,length(t)); % Allocate integrated W vector
for i = 1:length(t)-1
    W(i+1) = W(i)+sqrt(dt)*exp(th*t(i))*randn; 
end
ex = exp(-th*t);
x = starting_point*ex+target_point*(1-ex)+sig*ex.*W;
% figure;
% plot(t,x);