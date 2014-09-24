function plan = bipedGaitedMISOCP(robot, seed_plan, weights, goal_pos)
% @param seed_plan a blank footstep plan, provinding the structure of the
%                  desired plan. Probably generated with
%                  FootstepPlan.blank_plan()
% @param weights not used, maintained for compatibility with other planners
% @param goal_pos a struct with fields 'right' and 'left'.
%                 goal_pos.right is the desired 6 DOF pose
%                 of the right foot sole, and likewise for
%                 goal_pos.left

checkDependency('yalmip');

seed_plan.sanity_check();
rangecheck(seed_plan.footsteps(1).pos(6), -pi, pi);
rangecheck(seed_plan.footsteps(2).pos(6), -pi, pi);

% Temporary
assert(seed_plan.footsteps(1).frame_id == robot.foot_frame_id.right);

nsteps = length(seed_plan.footsteps);

seed_steps = [seed_plan.footsteps.pos];
prob = GaitedFootstepPlanningProblem();
prob.start_pose = struct('right', seed_steps(:,1), 'body', mean(seed_steps(:,1:2), 2), 'left', seed_steps(:,2));
prob.goal_pose = struct('right', goal_pos.right, 'left', goal_pos.left);
prob = prob.addIRISRegions(seed_plan.safe_regions);
prob.nframes = nsteps;
sol = prob.solveYalmip();
plan = sol.getBipedFootstepPlan(robot, seed_plan);

end