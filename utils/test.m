clear
clc
close all

fprintf("== Test for orthogonality == \n")
rng(2021)
matrix1 = eye(3) + 0.001*rand(3, 3);
matrix2 = Orientations.correctOrthogonality(matrix1);
fprintf("Initial matrix: \n")
disp(matrix1)
fprintf("Corrected matrix: \n")
disp(matrix2)

fprintf("== Test for Euler angles == \n")
eulerAngles1 = deg2rad([10; 10; 20]);
matrix1 = Orientations.eulerAnglesToRotation(eulerAngles1);
eulerAngles2 = Orientations.rotationToEulerAngles(matrix1);
fprintf("Euler angles [deg]: \n")
disp(rad2deg(eulerAngles1.'))
fprintf("Rotation matrix: \n")
disp(matrix1)
fprintf("Euler angles converted from the rotation matrix [deg]: \n")
disp(rad2deg(eulerAngles2.'))

fprintf("== Test for axis-angle representation == \n")
[a, phi] = Orientations.rotationToAxisAngle(matrix1);
matrix2 = Orientations.axisAngleToRotation(a, phi);
fprintf("Rotation matrix: \n")
disp(matrix1)
fprintf("Axis of rotation: \n")
disp(a.')
fprintf("Angle of rotation [deg]: \n")
disp(rad2deg(phi))
fprintf("Rotation matrix converted from the axis-angle: \n")
disp(matrix2);

fprintf("== Test for quaternion representation == \n")
a = [1/2; 1/2; 1/sqrt(2)];
phi = pi/3;
q = Orientations.axisAngleToQuat(a, phi);
R_quat = Orientations.quatToRotation(q);
R_axis = Orientations.axisAngleToRotation(a, phi);
q_re = Orientations.rotationToQuat(R_quat);
[a_re, phi_re] = Orientations.quatToAxisAngle(q_re);

fprintf("The axis of rotation: \n")
disp(a.')
fprintf("The angle of rotation: \n")
disp(phi)
fprintf("Quaternion: \n")
disp(q.')
fprintf("Rotation matrix converted from the quaternion: \n")
disp(R_quat)
fprintf("Rotation matrix converted from the axis-angle representation: \n")
disp(R_axis)
fprintf("Quaternion converted from the rotation matrix: \n")
disp(q_re.')
fprintf("The axis of rotation converted from the quaternion: \n")
disp(a_re.')
fprintf("The angle of rotation converted from the quaternion: \n")
disp(phi_re)