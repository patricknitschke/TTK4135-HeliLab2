%% Problem 2
%{
% Plotting of open loop estimates

figure(1)
subplot(511)
stairs(t,u),grid
ylabel('u = p_c')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('Travel [rad]')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('Travel rate [rad]')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('Pitch [rad]')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('time [s]'),ylabel('Pitch rate [s]')

% Plotting of helicopter data
states_struct = load('states_q01.mat');
states = states_struct.ans;
time = states(1,:).';
p_c = states(2,:).';
travel = states(3,:).';
travel_rate_r = states(4,:).';
p = states(5,:).';
p_dot = states(6,:).';

figure(2)
subplot(511)
stairs(time,p_c,'LineWidth',1),grid
ylabel('u = p_c [rad]')
subplot(512)
plot(time,travel,'m','LineWidth',1),grid
ylabel('Travel [rad]')
subplot(513)
plot(time,travel_rate_r,'m','LineWidth',1),grid
ylabel('Travel rate [rad]')
subplot(514)
plot(time,p,'m','LineWidth',1),grid
ylabel('Pitch [rad]')
subplot(515)
plot(time,p_dot,'m','LineWidth',1),grid
xlabel('time [s]'),ylabel('Pitch rate [rad]')

%% Problem 3
%Plotting of helicopter data
states_struct = load('states_LQR1.mat');
states = states_struct.ans;
time = states(1,:).';
p_c = states(2,:).';
travel = states(3,:).';
travel_rate_r = states(4,:).';
p = states(5,:).';
p_dot = states(6,:).';

for i = 1:length(time)
    travel(i) = travel(i) + x1_0;
end

figure(3)
subplot(511)
stairs(time,p_c,'LineWidth',1),grid
ylabel('u = p_c')
subplot(512)
plot(time,travel,'m','LineWidth',1),grid
ylabel('Travel [rad]')
subplot(513)
plot(time,travel_rate_r,'m','LineWidth',1),grid
ylabel('Travel rate [rad]')
subplot(514)
plot(time,p,'m','LineWidth',1),grid
ylabel('Pitch [rad]')
subplot(515)
plot(time,p_dot,'m','LineWidth',1),grid
xlabel('time [s]'),ylabel('Pitch rate [rad]')

%}
%% Problem 4

% Plotting of open loop estimates
t = 0:delta_t:delta_t*(length(u1)-1);
figure(4)
subplot(5,2,[1,2])
stairs(t,u1),grid
ylabel('u_1 = p_c')
subplot(5,2,[3,4])
stairs(t,u2),grid
ylabel('u_2 = e_c')
subplot(5,2,5)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('Travel [rad]')
subplot(5,2,6)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('Travel rate [rad]')
subplot(5,2,7)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('Pitch [rad]')
subplot(5,2,8)
plot(t,x4,'m',t,x4','mo'),grid
ylabel('Pitch rate [rad]')
subplot(5,2,9)
plot(t,x5,'m',t,x5','mo'),grid
ylabel('Elevation [rad]')
subplot(5,2,10)
plot(t,x6,'m',t,x6','mo'),grid
ylabel('Elev. rate [rad]')
xlabel('time [s]')


% Plotting of helicopter data
states_struct = load('states_withFeedback.mat');
states = states_struct.ans;
time = states(1,:).';
p_c = states(2,:).';
e_c = states(3,:).';
travel = states(4,:).';
travel_rate_r = states(5,:).';
p = states(6,:).';
p_dot = states(7,:).';
e = states(8,:).';
e_dot = states(9,:).';

for i = 1:length(time)
    travel(i) = travel(i) + x1_0;
end

figure(5)
subplot(5,2,[1,2])
stairs(time(1:9999),p_c(1:9999),'LineWidth',1),grid
ylabel('u_1 = p_c')
subplot(5,2,[3,4])
stairs(time(1:9999),e_c(1:9999),'LineWidth',1),grid
ylabel('u_2 = e_c')
subplot(5,2,5)
plot(time(1:9999),travel(1:9999),'m','LineWidth',1.2),grid
ylabel('Travel [rad]')
subplot(5,2,6)
plot(time(1:9999),travel_rate_r(1:9999),'m','LineWidth',1),grid
ylabel('Travel rate [rad]')
subplot(5,2,7)
plot(time(1:9999),p(1:9999),'m','LineWidth',1.2),grid
ylabel('Pitch [rad]')
subplot(5,2,8)
plot(time(1:9999),p_dot(1:9999),'m','LineWidth',1),grid
ylabel('Pitch rate [rad]')
subplot(5,2,9)
plot(time(1:9999),e(1:9999),'m','LineWidth',1.2),grid
ylabel('Elevation [rad]')
subplot(5,2,10)
plot(time(1:9999),e_dot(1:9999),'m','LineWidth',1),grid
ylabel('Elev. rate [rad]')
xlabel('time [s]')

% Plotting the constraint
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;

n = 1000;
lambdas = linspace(0, 3.5, 1000);
c_opt = zeros(1000,1);
for i = 1:1000
    c_opt(i) = alpha * exp(-beta*(lambdas(i)-(lambda_t)).^2);
    % - alpha * exp(-beta*(lambdas(i)-(lambda_t-1.5)).^2) ;
end

figure(6)
plot(lambdas, c_opt, 'LineWidth', 1.2)
grid on
title('Elevation constraint versus travel')
xlabel('Travel [rad]'), ylabel('Elevation constraint [rad]')


