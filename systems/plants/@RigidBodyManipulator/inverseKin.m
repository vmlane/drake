function q = inverseKin(obj,q0,varargin)
%
% inverseKin(obj,q0,body1,pos1,body2,pos2,...,options)
%
% attempts to solve the optimization problem
%   min_q (q-q0)'*(q-q0)
%   subject to 
%       body1 is at pos1
%       body2 is at pos2
%       ....
%   using q0 as the initial guess
% 
% @param q0 the initial pos
% @param body1...bodyN the index to the body that should be constrained.
% *You may use zero to indicate the center of mass*.
% @param pos1...posN the desired position of the [0;0;0] point in the corresponding body
%   is pos is 6 elements, then it also constrains the orientation 
%   (note: it does not make sense to specify an orientation for the COM).
% @option q_nom  replaces the cost function with (q-q_nom)'*(q-q_nom).
% @default q_nom = q0
% @option Q  puts a weight on the cost function qtilde'*Q*qtilde
%   

% todo: support orientation
% todo: if mex is enabled, I could move this entire loop (includins snopt)
% down into C.  It would get rid of a lot of overhead.

% todo: support state constraints, even contact constraints, but put an
% option in to enable/disable (presumably those will be harder to push into mex)

if isstruct(varargin{end}) 
  options = varargin{end};
  varargin=varargin(1:end-1);
else
  options = struct();
end

% note: keeping typecheck/sizecheck to a minimum because this might have to
% run inside a dynamical system (so should be fast)

if isfield(options,'q_nom') q_nom = options.q_nom; else q_nom = q0; end
if isfield(options,'Q') Q = options.Q; else Q = eye(obj.num_q); end

N = length(varargin);
nF = 3*N+1;

% support input as bodies instead of body inds
for i=1:2:length(varargin)
  if (isa(varargin{i},'RigidBody')) varargin{i} = find(obj.body==varargin{i},1); end
end

% mex is not implemented for orientations yet
use_mex = true;
for i=2:2:length(varargin)  
  if length(varargin{i})==6
    use_mex = false;
  end
end

function [f,G] = ik(q)
  f = zeros(nF,1); G = zeros(nF,obj.num_q);
  f(1) = (q-q_nom)'*Q*(q-q_nom);
  G(1,:) = 2*(q-q_nom)'*Q;
  kinsol = doKinematics(obj,q,false,use_mex);
  i=1;j=2;
  while i<length(varargin)
    if (varargin{i}==0)
      [x,J] = getCOM(obj,q);
    else
      [x,J] = forwardKin(obj,kinsol,varargin{i},[0;0;0]);
    end
    f([j:j+2]) = x - varargin{i+1}(1:3);
    G([j:j+2],:) = J;
    j=j+3;
    
    if length(varargin{i+1})==6
      [x,J] = forwardKinRot(obj,kinsol,varargin{i});
      f([j:j+2]) = x - varargin{i+1}(4:6);
      G([j:j+2],:) = J;
      j=j+3;
    end
    i=i+2; 
  end
%  if isfield(options,'visualizer') options.visualizer.draw(0,[q;0*q]); drawnow; end
end

global SNOPT_USERFUN;
SNOPT_USERFUN = @ik;

[q,F,info] = snopt(q0,obj.joint_limit_min,obj.joint_limit_max,zeros(nF,1),[inf;zeros(nF-1,1)],'snoptUserfun');%,0,1,[],[],[],iGfun,jGvar);

if (info~=1) 
  [str,cat] = snoptInfo(info);
  warning('SNOPT:InfoNotOne',['SNOPT exited w/ info = ',num2str(info),'.\n',cat,': ',str,'\n  Check p19 of Gill06 for more information.']);  
end

end