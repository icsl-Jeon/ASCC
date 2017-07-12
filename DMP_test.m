
[w,Mu_s,Sigma_s]=DMP_learning(obs_path');

[x_repro0,xdot_repro0]=DMP_repro(w,Mu_s,Sigma_s,1,obs_path(1,:)',obs_path(1,:)',obs_path(end,:)');
plot3(x_repro0(1,1:end),x_repro0(2,1:end),x_repro0(3,1:end))

hold on

start_step=10;
[x_repro,xdot_repro]=DMP_repro(w,Mu_s,Sigma_s,start_step,obs_path(10,:)'+[2 2 2]',obs_path(1,:)',obs_path(end,:)');
plot3(x_repro(1,start_step:end),x_repro(2,start_step:end),x_repro(3,start_step:end))

disp('change')
