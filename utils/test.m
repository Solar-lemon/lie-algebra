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
