%% Implement LQR with Optimal Trajectory
template_problem_2;

Q_lqr = diag([20,5,1,2]);
R_lqr = diag(0.0001);

A_lqr = A1;
B_lqr = B1;
N_lqr = 0;

[K_lqr, S, eigenv] = dlqr(A_lqr,B_lqr,Q_lqr,R_lqr,N_lqr);

x_opt = timeseries([x1 x2 x3 x4],t);
u_opt = timeseries(u, t);
