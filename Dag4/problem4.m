% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
problem3; % Change this to the init file corresponding to your helicopter

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A_c4 = zeros(6);
A_c4(1,2) = 1;
A_c4(2,3) = -K_2;
A_c4(3,4) = 1;
A_c4(4,3) = -K_1*K_pp;
A_c4(4,4) = -K_1*K_pd;
A_c4(5,6) = 1;
A_c4(6,5) = -K_3*K_ep;
A_c4(6,6) = -K_3*K_ed;

B_c4 = zeros(6,2);
B_c4(4,1) = K_1*K_pp;
B_c4(6,2) = K_3*K_ep;

A14 = eye(6) + delta_t * A_c4;
B14 = delta_t * B_c4;

% Number of states and inputs
mx = size(A14,2); % Number of states (number of columns in A)
mu = size(B14,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                              % p_dot
x5_0 = 0;                              % e
x6_0 = 0;                              % e_dot

x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';           % Initial values

% Time horizon and initialization
N  = 40;                                  % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;


% Bounds
ul      = -Inf*ones(mu,1);              % Lower bound on states (no bound)
uu      = Inf*ones(mu,1);               % Upper bound on states (no bound)
ul(1)   = -30*pi/180;                   % Lower bound on control
uu(1) 	= 30*pi/180;                   % Upper bound on control
%e bounds
xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on state x3
xu(3)   = uu(1);                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;
Q1(5,5) = 0;
Q1(6,6) = 0;

% Weight on state x4
P1 = zeros(mu,mu);
P1(1,1) = 0.5;                                 % Weight on input
P1(2,2) = 0.5;

Q4 = gen_q(Q1,P1,N,M);                   % Generate Q, hint: gen_q
c4 = zeros(N*mx+M*mu,1);                 % Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model
Aeq = gen_aeq(A14,B14,N,mx,mu);         % Generate A, hint: gen_aeq

A0 = A14;
beq_1 = A0*x0;
beq = Aeq*z0; % Generate b


beq(1) = beq_1(1);
z0(1) = x0(1);                                 % Initial value for optimization


%% Solve optimisation problem with nonlinear constraints
tic
%[z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub,x0); % hint: quadprog. Type 'doc quadprog' for more info

options = optimoptions('fmincon','Algorithm','sqp');
min_func = @(z) z'*Q4*z;

[z,lambda] = fmincon(min_func,z0,[],[], Aeq, beq, vlb, vub, @nonlcon_e,options);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q4(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u1  = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu-1)]; % Control input from solution
u2  = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution


num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2  = [zero_padding; u2; zero_padding];
x1  = [x1_0*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

%% Plotting of open loop estimates
t = 0:delta_t:delta_t*(length(u1)-1);

figure(8)
subplot(811)
stairs(t,u1),grid
ylabel('p_c')
subplot(812)
stairs(t,u2),grid
ylabel('e_c')
subplot(813)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','mo'),grid
ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5','mo'),grid
ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6','mo'),grid
ylabel('edot')
xlabel('tid (s)')


%% Plotting of helicopter data

states_struct = load('states.mat');
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

figure(9)
subplot(811)
stairs(time,p_c,'LineWidth',1.2),grid
ylabel('p_c')
subplot(812)
stairs(time,e_c,'LineWidth',1.2),grid
ylabel('e_c')
subplot(813)
plot(time,travel,'m','LineWidth',1.2),grid
ylabel('lambda')
subplot(814)
plot(time,travel_rate_r,'m','LineWidth',1.2),grid
ylabel('r')
subplot(815)
plot(time,p,'m','LineWidth',1.2),grid
ylabel('p')
subplot(816)
plot(time,p_dot,'m','LineWidth',1.2),grid
ylabel('pdot')
subplot(817)
plot(time,e,'m','LineWidth',1.2),grid
ylabel('e')
subplot(818)
plot(time,e_dot,'m','LineWidth',1.2),grid
ylabel('edot')
xlabel('tid (s)')



%% LQR
Q_lqr4 = diag([20,5,1,2,1,1]);
R_lqr4 = diag([0.1 0.1]);

A_lqr4 = A14;
B_lqr4 = B14;
N_lqr4 = 0;

[K_lqr, S4, eigenv4] = dlqr(A_lqr4,B_lqr4,Q_lqr4,R_lqr4,N_lqr4);

x_opt = timeseries([x1 x2 x3 x4 x5 x6], t);
u_opt = timeseries([u1 u2], t);









