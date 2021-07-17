classdef Orientations < handle
    properties(Constant)
        X_AXIS = 1;
        Y_AXIS = 2;
        Z_AXIS = 3;
    end
    methods(Static)
        function isOrthogonal = checkOrthogonality(matrix)
            isOrthogonal = (norm(matrix*matrix.' - eye(3)) < 1e-4);
        end
        
        function matrix = correctOrthogonality(matrix)
            matrix = real(sqrtm(matrix*matrix.')) \ matrix;
        end
        
        function isSkewSymmetric = checkSkewSymmetry(matrix)
            isSkewSymmetric = (norm(matrix + matrix.') < 1e-4);
        end
        
        function matrix = correctSkewSymmetry(matrix)
            matrix = (matrix - matrix.')/2;
        end
        
        function isUnity = checkUnity(v)
            isUnity = (abs(norm(v) - 1) < 1e-4);
        end
        
        function v = correctUnity(v)
            v = v/norm(v);
        end
        
        function matrix = basicRotation(axis, phi)
            % Suppose that the vehicle frame {v} is rotated from
            % the inertial frame {i} by theta with respect to X_AXIS.
            % Then the rotation matrix R_{vi} from {i} to {v} can be
            % expressed as
            % R_{vi} = [1, 0, 0;
            %           0, cos(theta), sin(theta);
            %           0, -sin(theta), cos(theta)];
            % and pos_{v} = R_{vi}*pos_{i}
            switch axis
                case Orientations.X_AXIS
                    matrix = [...
                        1, 0, 0;
                        0, cos(phi), sin(phi);
                        0, -sin(phi), cos(phi)];
                case Orientations.Y_AXIS
                    matrix = [...
                        cos(phi), 0, -sin(phi);
                        0, 1, 0;
                        sin(phi), 0, cos(phi)];
                case Orientations.Z_AXIS
                    matrix = [...
                        cos(phi), sin(phi), 0;
                        -sin(phi), cos(phi), 0;
                        0, 0, 1];
            end
        end
        
        function matrix = hat(v)
            % v: 3x1 vector
            matrix = [...
                0, -v(3), v(2);
                v(3), 0, -v(1);
                -v(2), v(1), 0];
        end
        
        function matrix = axisAngleToRotation(a, phi)
            % a: 3x1 unit vector, axis of rotation
            % phi: scalar, angle of rotation
            % the resulting rotation matrix is from inertial frame
            % to body frame.
            a_hat = Orientations.hat(a);
            matrix = eye(3) - sin(phi)*a_hat + (1 - cos(phi))*a_hat^2;
        end
        
        function [a, phi] = rotationToAxisAngle(matrix)
            % phi: angle of rotation from inertial frame to body frame
            % a: axis of rotation from inertial frame to body frame
            [eigenVectors, eigenValues] = eig(matrix);
            differences = abs(diag(eigenValues) - 1);
            [minDifference, index] = min(differences);
            assert(minDifference < 1e-4, "No eigenvalue exists that is close to 1.")
            a = eigenVectors(:, index);
            a = real(a);
            phi = real(acos((trace(matrix) - 1)/2));
            
            if norm(Orientations.axisAngleToRotation(a, phi) - matrix) > 1e-4
                a = -a;
                assert(...
                    norm(Orientations.axisAngleToRotation(a, phi) - matrix) < 1e-4,...
                    "No solution exists")
            end
        end
        
        function matrix = eulerAnglesToRotation(varargin)
            % convention: 3-2-1 Euler angle
            % (the resulting rotation matrix is from inertial frame
            % to body frame)
            % varargin = {eulerAngles} where eulerAngles = [phi; theta; psi] or
            % varargin = {phi, theta, psi}
            switch numel(varargin)
                case 1
                    eulerAngles = varargin{1};
                    phi   = eulerAngles(1);
                    theta = eulerAngles(2);
                    psi   = eulerAngles(3);
                case 3
                    phi   = varargin{1};
                    theta = varargin{2};
                    psi   = varargin{3};
            end
            R_psi = [...
                cos(psi), sin(psi), 0;
                -sin(psi), cos(psi), 0;
                0, 0, 1];
            R_theta = [...
                cos(theta), 0, -sin(theta);
                0, 1, 0;
                sin(theta), 0, cos(theta)];
            R_phi = [...
                1, 0, 0;
                0, cos(phi), sin(phi);
                0, -sin(phi), cos(phi)];
            matrix = R_phi*R_theta*R_psi;
        end
        
        function eulerAngles = rotationToEulerAngles(matrix)
            r_11 = matrix(1, 1);
            r_12 = matrix(1, 2);
            r_13 = matrix(1, 3);
            r_23 = matrix(2, 3);
            r_33 = matrix(3, 3);
            
            if abs(r_13) < 1 - 1e-6
                theta = atan2(-r_13, sqrt(r_11^2 + r_12^2));
                psi = atan2(r_12, r_11);
                phi = atan2(r_23, r_33);
            else
                r_21 = obj.matrix(2, 1);
                r_22 = obj.matrix(2, 2);
                if r_13 < 0
                    theta = pi/2;
                    psi = 0;
                    phi = atan2(r_21, r_22);
                else
                    theta = -pi/2;
                    psi = 0;
                    phi = -atan2(r_21, r_22);
                end
            end
            eulerAngles = [phi; theta; psi];
        end
        
        function q = axisAngleToQuat(a, phi)
            q = [cos(phi/2); sin(phi/2)*a];
        end
        
        function [a, phi] = quatToAxisAngle(q)
            a = q(2:4);
            if norm(a) > 1e-8
                a = a/norm(a);
            end
            q_1 = max(-1, min(1, q(1)));
            phi = 2*acos(q_1);
        end
        
        function R = quatToRotation(q)
            % q: 4x1 (elements of a quaternion)
            eta = q(1);
            epsilon = q(2:4);
            R = (eta^2 - epsilon.'*epsilon)*eye(3) + ...
                2*epsilon*(epsilon.') - 2*eta*Orientations.hat(epsilon);
        end
        
        function q = rotationToQuat(R)
            % R: 3x3 matrix
            r_11 = R(1, 1);
            r_22 = R(2, 2);
            r_33 = R(3, 3);
            
            r_12 = R(1, 2);
            r_21 = R(2, 1);
            r_23 = R(2, 3);
            r_32 = R(3, 2);
            r_31 = R(3, 1);
            r_13 = R(1, 3);
            
            eta = 1/2*sqrt(...
                max(0, 1 + trace(R)));
            epsilon_1 = sign(r_23 - r_32)*1/2*sqrt(...
                max(0, 1 + r_11 - r_22 - r_33));
            epsilon_2 = sign(r_31 - r_13)*1/2*sqrt(...
                max(0, 1 - r_11 + r_22 - r_33));
            epsilon_3 = sign(r_12 - r_21)*1/2*sqrt(...
                max(0, 1 - r_11 - r_22 + r_33));
            
            q = [eta; epsilon_1; epsilon_2; epsilon_3];
        end
    end
end
