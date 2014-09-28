function testQuadrupedPlanner()
% Use the gaited mixed-integer footstep planner to solve a quadruped planning problem. 

prob = GaitedFootstepPlanningProblem();

% Pre-generate some safe regions in IRIS region format
safe_regions = struct('A', {}, 'b', {}, 'point', {}, 'normal', {});
V = [-0.15, 0.2 .2, -0.15; 0.2, 0.2, -0.2, -0.2];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.3;0;1], 'normal', [0;0;1]);

V = [0.25, 0.4 .4, 0.25; 0, 0, -0.3, -0.3];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.3;0;1], 'normal', [0;0;1]);

V = [.59, 2, 2, .59; 1, 1, -1, -1];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.8;0;1], 'normal', [0;0;1]);

V = [.5,.51,.5,.51; -.05, -.05, -.06, -.06];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.5;-.05;1], 'normal', [0;0;1]);

V = [.5,.51,.51,.5; .05, .05, .06, .06];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.5;.05;1], 'normal', [0;0;1]);

V = [.45,.46,.46,.45; 0, 0, .01, .01];
[A, b] = poly2lincon(V(1,:), V(2,:));
safe_regions(end+1) = struct('A', A, 'b', b, 'point', [.5;.05;1], 'normal', [0;0;1]);

% Convert safe regions from x,y to x,y,yaw
for j = 1:length(safe_regions)
  safe_regions(j).A = [safe_regions(j).A, zeros(size(safe_regions(j).A, 1), 1)];
end

prob = prob.addIRISRegions(safe_regions);

prob.feet = {'rf', 'lf', 'rh', 'lh'};
prob.foci = struct('rf', struct('v', {[0.1; -0.05]}, 'r', {0.05}),...
              'lf', struct('v', {[0.1; 0.05]}, 'r', {0.05}),...
              'rh', struct('v', {[-0.1; -0.05]}, 'r', {0.05}),...
              'lh', struct('v', {[-0.1; 0.05]}, 'r', {0.05}));

% Limits on delta Z and yaw from body to each foot
A = [0,0,1,0,0,0;
     0,0,-1,0,0,0;
     0,0,0,0,0,1;
     0,0,0,0,0,-1];
z_range = [-.3, -.1];
yaw_range = [-pi, pi];
b = [z_range(2);-z_range(1);yaw_range(2);-yaw_range(1)];
% A = [0,0,1,0,0,0;
%      0,0,-1,0,0,0];
% b = [.2;.2];
lcon_struct = struct('A', A, 'b', b);
prob.body_to_feet_constraints = struct('rf', lcon_struct,...
                                       'lf', lcon_struct,...
                                       'rh', lcon_struct,...
                                       'lh', lcon_struct);
prob.gait = struct('rf', {}, 'lf', {}, 'rh', {}, 'lh', {});
prob.gait(1) = struct('rf', true, 'lf', true, 'rh', true, 'lh', true);
prob.gait(2) = struct('rf', true, 'lf', true, 'rh', false, 'lh', true);
prob.gait(3) = struct('rf', true, 'lf', false, 'rh', false, 'lh', true);
prob.gait(4) = struct('rf', true, 'lf', false, 'rh', true, 'lh', true);
prob.gait(5) = struct('rf', true, 'lf', true, 'rh', true, 'lh', true);
prob.gait(6) = struct('rf', true, 'lf', true, 'rh', true, 'lh', false);
prob.gait(7) = struct('rf', false, 'lf', true, 'rh', true, 'lh', false);
prob.gait(8) = struct('rf', false, 'lf', true, 'rh', true, 'lh', true);

prob.nframes = 32;
prob.start_pose = struct('body', [0;0;1.2;0;0;0],...
                         'rf', [0.1;-0.05;0;0;0;0],...
                         'lf', [0.1;0.05;0;0;0;0],...
                         'rh', [-0.1;-0.05;0;0;0;0],...
                         'lh', [-0.1;0.05;0;0;0;0]);
prob.goal_pose = struct('body', [1;0;0.2;0;0;0]);

sol = prob.solveYalmip();
% save('sol.mat', 'sol');
[xtraj, v] = sol.getSimpleGaitTrajectory();
v.playback(xtraj, struct('slider', true));