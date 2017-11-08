function spline_sample = sample_spline(spline_variance,obs_variance,ctrl_pts,obs)


t = linspace(0,1,size(obs,1))';
spline.cov = spline_variance .* eye(length(t).*2);
obs_cov = obs_variance .* eye(length(t).*2);

a = 1;
subA = [(1-a.*t).^3,3.*(1-a.*t).^2.*a.*t,3.*(1-a.*t).*a.*t.^2,a.*t.^3];
A = blkdiag(subA,subA);


obs_cov = obs_variance .* eye(length(t).*2);
distribution.traj.cov = inv((inv(spline.cov) + inv(obs_cov)));
distribution.traj.mean = distribution.traj.cov * (inv(spline.cov) * A * reshape(ctrl_pts,[],1) + inv(obs_cov) * reshape(obs,[],1));
distribution.traj.mean = reshape(distribution.traj.mean,[],2);

spline_sample = reshape(mvnrnd(reshape(distribution.traj.mean,[],1),distribution.traj.cov),[],2);