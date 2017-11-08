function distribution = learn_spline(obs)

t = linspace(0,1,size(obs,1))';

a = 1;
variance = 1e2;
gamma = 1e0;

subA = [(1-a.*t).^3,3.*(1-a.*t).^2.*a.*t,3.*(1-a.*t).*a.*t.^2,a.*t.^3];
A = blkdiag(subA,subA);

spline.cov = gamma .* eye(length(t).*2);
distribution.ctrl_pts.cov = inv(A' * inv(spline.cov) * A);
distribution.ctrl_pts.mean = reshape(distribution.ctrl_pts.cov * A' * inv(spline.cov) * reshape(obs,[],1),[],2);

obs_cov = variance .* eye(length(t).*2);
distribution.traj.cov = inv((inv(spline.cov) + inv(obs_cov)));
distribution.traj.mean = distribution.traj.cov * (inv(spline.cov) * A * reshape(distribution.ctrl_pts.mean,[],1) + inv(obs_cov) * reshape(obs,[],1));
distribution.traj.mean = reshape(distribution.traj.mean,[],2);

