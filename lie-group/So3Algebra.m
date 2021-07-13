classdef So3Algebra < handle
    properties
        vector
        matrix
    end
    methods
        function obj = So3Algebra(inValue)
            if all(size(inValue) == [3, 1])
                initializeWithVector(obj, inValue);
                return
            end
            if all(size(inValue) == [3, 3])
                initializeWithMatrix(obj, inValue);
                return
            end
        end
        
        function [a, phi] = axisAngle(obj)
            phi = norm(obj.vector);
            a = obj.vector/phi;
        end
        
        function out = bracket(obj, otherObj)
            newVector = obj.matrix*otherObj.vector;
            out = So3Algebra(newVector);
        end
        
        function out = exp(obj)
            phi = norm(obj.vector);
            if phi < 1e-8
                out = SO3Group(eye(3));
                return
            end
            a_hat = obj.matrix/phi;
            expMatrix = eye(3) + sin(phi)*a_hat + (1 - cos(phi))*a_hat^2;
            out = SO3Group(expMatrix);
        end
        
        function out = jacobian(obj)
            out = leftJacobian(obj);
        end

        function out = leftJacobian(obj)
            phi = norm(obj.vector);
            if phi < 1e-4
                out = eye(3);
                return
            end
            a_hat = obj.matrix/phi;
            out = eye(3) + (1 - cos(phi))/phi*a_hat + (1 - sin(phi)/phi)*a_hat^2;
        end
        
        function out = rightJacobian(obj)
            phi = norm(obj.vector);
            if phi < 1e-4
                out = eye(3);
                return
            end
            a_hat = obj.matrix/phi;
            out = eye(3) - (1 - cos(phi))/phi*a_hat + (1 - sin(phi)/phi)*a_hat^2;
        end
        
        function out = approximateLeftJacobian(obj, N)
            temp = eye(3);
            out = temp;
            for i = 1:N
                temp = 1/(i + 1)*obj.matrix*temp;
                out = out + temp;
            end
        end
        
        function out = approximateRightJacobian(obj, N)
            temp = eye(3);
            out = temp;
            for i = 1:N
                temp = -1/(i + 1)*obj.matrix*temp;
                out = out + temp;
            end
        end
        
        function out = inverseJacobian(obj)
            phi = norm(obj.vector);
            if phi < 1e-4
                out = eye(3);
                return
            end
            a_hat = obj.matrix/phi;
            out = eye(3) - phi/2*a_hat + (1 - phi/2*cot(phi/2))*a_hat^2;
        end
        
        function out = detJacobian(obj)
            phi = norm(obj.vector);
            if phi < 1e-4
                out = 1;
                return
            end
            out = 2*(1 - cos(phi))/phi^2;
        end
    end
    
    methods(Access = private)
        function initializeWithVector(obj, inValue)
            obj.vector = inValue;
            obj.matrix = So3Algebra.vectorToMatrix(inValue);
        end
        
        function initializeWithMatrix(obj, inValue)
            obj.vector = So3Algebra.matrixToVector(inValue);
            obj.matrix = inValue;
        end
    end
    
    methods(Static)
        function matrix = vectorToMatrix(vector)
            matrix = [...
                0, -vector(3), vector(2);
                vector(3), 0, -vector(1);
                -vector(2), vector(1), 0];
        end
        
        function vector = matrixToVector(matrix)
            isSkewSymmetric = Orientations.checkSkewSymmetry(matrix);
            assert(isSkewSymmetric, "A matrix for So3Algebra has to be skew-symmetric.")
            
            matrix = Orientations.correctSkewSymmetry(matrix);
            vector = [matrix(3, 2); -matrix(3, 1); matrix(2, 1)];
        end
    end
end