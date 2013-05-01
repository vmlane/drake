function testOrientation
% test rpy to rotmat
for i = 1:1000
    rpy = 2*pi*(rand(3,1)-0.5);
    rpy1 = rotmat2rpy(rpy2rotmat(rpy));
    mat = rpy2rotmat(rpy);
    mat1 = rpy2rotmat(rpy1);
    valuecheck(mat,mat1,1e-6);
end
disp('rpy to rotmat and rotmat to rpy is consistent');
% test quat to rotmat
for i = 1:1000
    quat = randn(4,1);
    quat = quat./norm(quat);
    mat = quat2rotmat(quat);
    quat1 = rotmat2quat(mat);
    mat1 = quat2rotmat(quat1);
    valuecheck(mat,mat1,1e-6);
end
disp('quat to rotmat and rotmat to quat is consistent');

for i = 1:1000
    rpy = 2*pi*(rand(3,1)-0.5);
    quat = rpy2quat(rpy);
    quat1 = rotmat2quat(rpy2rotmat(rpy));
    valuecheck(1-abs(quat1'*quat),0,1e-6);
    rpy1 = quat2rpy(quat);
    mat = rpy2rotmat(rpy);
    mat1 = rpy2rotmat(rpy1);
    valuecheck(mat,mat1,1e-6);
end
disp('rpy to quat and quat to rpy is consistent, rpy to quat is correct');

for i = 1:1000
    quat = randn(4,1);
    quat = quat./norm(quat);
    axis = quat2axis(quat);
    quat1 = axis2quat(axis);
    valuecheck(acos(abs(quat'*quat1)),0,1e-6);
end
disp('axis to quat and quat to axis is consistent');

for i = 1:1000
    quat = randn(4,1);
    quat = quat./norm(quat);
    mat = quat2rotmat(quat);
    axis = rotmat2axis(mat);
    quat1 = axis2quat(axis);
    valuecheck(acos(abs(quat'*quat1)),0,1e-6);
end
disp('rotmat2axis is correct');
end