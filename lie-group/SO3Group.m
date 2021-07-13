classdef SO3Group < handle
    % Special orthogonal group SO(3) is the set of valid rotation matrices
    % For an element C in SO(3), the following property holds.
    % CC^T = 1, detC = 1
    % Generally, an element of SO(3) group represents a rotation matrix
    % from body frame to inertial frame.
    properties
        matrix
    end
    
    methods
        function obj = SO3Group(matrix)
            isOrthogonal = Orientations.checkOrthogonality(matrix);
            assert(isOrthogonal, "A rotation matrix has to be orthogonal.")
            
            obj.matrix = Orientations.correctOrthogonality(matrix);
        end
        
        function out = inverse(obj)
            out = SO3Group(obj.matrix.');
        end
        
        function out = log(obj)
            % Generally, phi represents the angle of rotation from inertial
            % frame to body frame and a represents the axis of rotation.
            [eigenVectors, eigenValues] = eig(obj.matrix);
            differences = abs(diag(eigenValues) - 1);
            [minDifference, index] = min(differences);
            assert(minDifference < 1e-4, "No eigenvalue exists that is close to 1.")
            a = eigenVectors(:, index);
            a = real(a);
            phi = real(acos((trace(obj.matrix) - 1)/2));
            
            vector = phi*a;
            if norm(So3Algebra(vector).exp().matrix - obj.matrix) > 1e-4
                a = -a;
                vector = phi*a;
                assert(norm(So3Algebra(vector).exp().matrix - obj.matrix) < 1e-4,...
                    "No solution exists")
            end
            out = So3Algebra(vector);
        end
        
        function out = difference(obj, otherObj)
            % Using left difference,
            % the difference of C_2 relative to C_1 can be defined as
            % epsilon_hat = ln(C_2*inv(C_1)) = ln(C_2*C_1^T)
            % which means C_2 = exp(epsilon_hat)*C_1
            out = SO3Group(otherObj.matrix * (obj.matrix.')).log();
        end
    end
    
    methods
        % operator overloading
        function out = mtimes(obj1, obj2)
            out = SO3Group(obj1.matrix*obj2.matrix);
        end
        
        function out = mldivide(obj1, obj2)
            out = obj1.matrix.'*obj2.matrix;
        end
        
        function out = mrdivide(obj1, obj2)
            out = obj1.matrix*(obj2.matrix.');
        end
    end
    
    methods(Static)
        function test()
            clc
            close all
            
            fprintf("== Test for SO3Group == \n")
            % matrix: a rotation matrix from inertial frame to body frame.
            matrix = Orientations.basicRotation(Orientations.X_AXIS, deg2rad(30));
            element = SO3Group(matrix.');
            [axis, phi] = element.log().axisAngle();
            
            fprintf("rotation matrix: \n")
            disp(element.matrix)
            fprintf("rotation axis: \n")
            disp(axis.')
            fprintf("rotation angle: %.1f[deg] \n", rad2deg(phi))
        end
    end
end