clear all; close all; clc

pos(:,1) = [5;-10]; 
pos(:,2) = [18;18]; 
pos(:,3) = [38;-5]; 
pos(:,4) = [45;15];
a = 1;
variance = 1e2;
gamma = 1e0;

color = [[117,112,179]./255; [217,95,2]./255; 0.7,0.7,0.7];

t = linspace(0,2,100)';

subA = [(1-a.*t).^3,3.*(1-a.*t).^2.*a.*t,3.*(1-a.*t).*a.*t.^2,a.*t.^3];
A = blkdiag(subA,subA);

trajectory.true = reshape(A * reshape(pos',[],1),[],2);
trajectory.obs =  trajectory.true + sqrt(variance) .* randn(size(trajectory.true));

figure(1)
plot(trajectory.true(:,1),trajectory.true(:,2),'LineWidth',2,'Color',color(2,:))
hold on; plot(trajectory.obs(:,1),trajectory.obs(:,2),'*','Color',color(2,:))

spline.cov = gamma .* eye(length(t).*2);
distribution.ctrl_pts.cov = inv(A' * inv(spline.cov) * A);
distribution.ctrl_pts.mean = reshape(distribution.ctrl_pts.cov * A' * inv(spline.cov) * reshape(trajectory.obs,[],1),[],2);

obs.cov = variance .* eye(length(t).*2);
distribution.traj.cov = inv((inv(spline.cov) + inv(obs.cov)));
distribution.traj.mean = distribution.traj.cov * (inv(spline.cov) * A * reshape(distribution.ctrl_pts.mean,[],1) + inv(obs.cov) * reshape(trajectory.obs,[],1));
distribution.traj.mean = reshape(distribution.traj.mean,[],2);

hold on; plot(distribution.traj.mean(:,1),distribution.traj.mean(:,2),'LineWidth',2,'Color',color(1,:))
legend('true trajectory','observations','learned trajectory')