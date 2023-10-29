classdef RwArray < handle

    properties(SetAccess = protected, GetAccess = public)
        reactionwheels           % array of reaction wheels(X, Y, and Z unless set up otherwise)
        installationMatrix       % each column corresponds to rw axis direction in the body frame
        torqueAllocationMetric   % matrix of the torque quadratic form to be minimized at torque allocation
        aInv                     % the iverse of the torqueAllocationMetric
    end

    methods
        function this = RwArray(parameters)
            arguments
                parameters.baselineRw ReactionWheel
                parameters.rwConfiguration {mustBeA(parameters.rwConfiguration, 'RwArrayConfiguration')} = RwArrayConfiguration.standardXYZ;
                parameters.dcmArray(3, 3, :) {mustBeNumeric} = []; % TODO
            end

            switch parameters.rwConfiguration
                case RwArrayConfiguration.standardXYZ
                    rwCount = 3;
                case RwArrayConfiguration.regularTetrahedron
                    rwCount = 4;
                case RwArrayConfiguration.userDefined
                    rwCount = size(parameters.dcmArray, 3);
            end

            rwArray(1, rwCount) = ReactionWheel();
            for rwIdx = 1:rwCount
                rwArray(rwIdx).setMaxTorque(parameters.baselineRw.maxTorque);
                rwArray(rwIdx).setMaxAngMomentum(parameters.baselineRw.maxAngMomentum);
            end

            this.reactionwheels = rwArray;

            this.installationMatrix = zeros(3, rwCount);
            this.setupRwConfiguration(parameters.rwConfiguration, parameters.baselineRw.axis)

            this.torqueAllocationMetric = eye(rwCount);
            this.aInv = inv(this.torqueAllocationMetric);
        end

        function actuatedTorque = actuateCommand(this, commandTorque, duration)
            actuatedTorque = zeros(3, 1);
            torqueAllocation = this.allocateTorque(commandTorque);

            for rwIdx = 1:size(this.reactionwheels, 2)
                actuatedTorque = actuatedTorque + ...
                                 this.reactionwheels(rwIdx).actuateControlTorque(torqueAllocation(rwIdx), ...
                                                                                 duration);
            end
        end
    end

    methods(Access = private)

        function h = angularMomentumInArray(this)
            h = [this.reactionwheels.angularMomentum]';
        end        

        function torqueAllocation = allocateTorque(this, commandTorque)
            % TODO do something better for allocation

            Ainv = this.aInv;
            WT = this.installationMatrix';
            W = this.installationMatrix;

            torqueAllocation = (Ainv * WT / (W * Ainv * WT)) *  commandTorque;
        end

        function setupRwConfiguration(this, configurationType, defaultAxis)
            switch configurationType
                case RwArrayConfiguration.standardXYZ
                    for rwAllocation = ['X', 'Y', 'Z']
                        this.update3rwConfiguration(rwAllocation, defaultAxis);
                    end
                case RwArrayConfiguration.regularTetrahedron
                    for rwAllocation = ['A', 'B', 'C', 'D']
                        this.updateTetrahedronConfiguration(rwAllocation, defaultAxis);
                    end
                case RwArrayConfiguration.userDefined
                    %TODO
            end
        end

        function update3rwConfiguration(this, rwAllocation, defaultAxis)
            switch rwAllocation
                case 'X'
                    rwIdx = 1;
                    rwDcm = [cos(pi / 2) 0 sin(pi / 2)
                                  0       1     0
                             -sin(pi / 2) 0 cos(pi / 2)];
                case 'Y'
                    rwIdx = 2;
                    rwDcm = [1       0             0
                              0   cos(-pi / 2)  -sin(-pi / 2)
                              0   sin(-pi / 2)   cos(-pi / 2)];
                case 'Z'
                    rwIdx = 3;
                    rwDcm = eye(3);
            end

            this.reactionwheels(rwIdx).setName(['rw' rwAllocation]);
            this.reactionwheels(rwIdx).setDcm(rwDcm);
            installatioAxis = rwDcm * defaultAxis;
            this.installationMatrix(:, rwIdx) = installatioAxis / vecnorm(installatioAxis);
        end

        function updateTetrahedronConfiguration(this, rwAllocation, defaultAxis)
            function rwDcm = calcDcm(vertexAxis, originalAxis)
                vecFrom  = originalAxis / vecnorm(originalAxis);
                vecTo = vertexAxis / vecnorm(vertexAxis);

                vecAngle = atan2(vecnorm(crossProduct(vecTo, vecFrom)), dotProduct(vecTo, vecFrom));
                quatAxis = crossProduct(vecTo, vecFrom);
                if vecnorm(quatAxis) > 1e-7
                    quatAxis = quatAxis  / vecnorm(quatAxis);
                end

                rwDcm = quaternion2dcm([cos(vecAngle / 2); quatAxis * sin(vecAngle / 2)]');
            end

            switch rwAllocation
                case 'A'
                    rwIdx = 1;
                    vertexA = [sqrt(8 / 9); 0; -1 / 3];
                    rwDcm = calcDcm(vertexA, defaultAxis);
                case 'B'
                    rwIdx = 2;
                    vertexB = [-sqrt(2 / 9); sqrt(2 / 3); -1 / 3];
                    rwDcm = calcDcm(vertexB, defaultAxis);
                case 'C'
                    rwIdx = 3;
                    vertexC = [-sqrt(2 / 9); -sqrt(2 / 3); -1 / 3];
                    rwDcm = calcDcm(vertexC, defaultAxis);
                case 'D'
                    rwIdx = 4;
                    vertexD = [0; 0; 1];
                    rwDcm = calcDcm(vertexD, defaultAxis);                    
            end

            this.reactionwheels(rwIdx).setName(['rw' rwAllocation]);
            this.reactionwheels(rwIdx).setDcm(rwDcm);
            installatioAxis = rwDcm * defaultAxis;
            this.installationMatrix(:, rwIdx) = installatioAxis / vecnorm(installatioAxis);
        end

    end
end