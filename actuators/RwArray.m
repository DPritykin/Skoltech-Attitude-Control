classdef RwArray < handle

    properties(SetAccess = protected, GetAccess = public)
        reactionwheels           % array of reaction wheels(X, Y, and Z unless set up otherwise)
    end

    methods
        function this = RwArray(parameters)
            arguments
                parameters.baselineRw ReactionWheel
                parameters.rwCount {mustBeNumericOrLogical} = 4;
                parameters.doStandardXyzArray {mustBeNumericOrLogical} = false;
            end

            rwArray(1, parameters.rwCount) = ReactionWheel();
            for rwIdx = 1:parameters.rwCount
                rwArray(rwIdx).setMaxTorque(parameters.baselineRw.maxTorque);
                rwArray(rwIdx).setMaxAngMomentum(parameters.baselineRw.maxAngMomentum);
            end

            this.reactionwheels = rwArray;

            if parameters.doStandardXyzArray
                for rwAxisBf = ['X', 'Y', 'Z']
                    this.updateReactionWheelStandard(rwAxisBf)
                end
            else
                for rwAxisBf = ['A', 'B', 'C', 'D']
                    this.updateReactionWheelPyramid(rwAxisBf)
                end
            end
        end


        function actuatedTorque = actuateCommand(this, commandTorque, duration, rwAngMomentum)
            
            WheelToBody = (1/sqrt(3))*[1 -1 -1 1; 1  1 -1 -1; 1 1 1 1]; 
            BodyToWheel = (sqrt(3)/4)*[1 1 1; -1 1 1; -1 -1 1; 1 -1 1];

            actuatedTorque = zeros(3, 1);
            
            if length(this.reactionwheels) == 4
                commandTorqueWheel = BodyToWheel * commandTorque;
                rwAngMomentumWheel = BodyToWheel * rwAngMomentum;
            else 
                commandTorqueWheel = commandTorque;
                rwAngMomentumWheel = rwAngMomentum;
            end 

            for rwIdx = 1:size(this.reactionwheels, 2)  % Iterate through the reaction wheels
                actuatedTorque = actuatedTorque + ...
                    this.reactionwheels(rwIdx).actuateControlTorque(commandTorqueWheel(rwIdx), duration, rwAngMomentumWheel(rwIdx));
            end
        end

    end

    methods(Access = private)

        function updateReactionWheelPyramid(this, rwAxisBf)
            switch rwAxisBf
                case 'A'
                    rwIdx = 1;
                    rwName = 'rwA';
                    rwDcm = eye(3);

                case 'B'
                    rwIdx = 2;
                    rwName = 'rwB';
                    rwDcm = [cos(2*pi/3), -sin(2*pi/3), 0;
                             sin(2*pi/3), cos(2*pi/3), 0;
                             0, 0, 1];
                  case 'C'  
                    rwIdx = 3;
                    rwName = 'rwC';
                    rwDcm = [cos(4*pi/3), -sin(4*pi/3), 0;
                             sin(4*pi/3), cos(4*pi/3), 0;
                             0, 0, 1];

                case 'D'
                    rwIdx = 4;
                    rwName = 'rwD';
                    rwDcm = [cos(pi), -sin(pi), 0;
                               sin(pi), cos(pi), 0;
                               0, 0, 1];
            end

            this.reactionwheels(rwIdx).setName(rwName);
            this.reactionwheels(rwIdx).setDcm(rwDcm);
        end

        function updateReactionWheelStandard(this, rwAxisBf)
            switch rwAxisBf
                case 'X'
                    rwIdx = 1;
                    rwName = 'rwX';
                    rwDcm = [cos(pi / 2) 0 sin(pi / 2)
                                  0       1     0
                             -sin(pi / 2) 0 cos(pi / 2)];
                case 'Y'
                    rwIdx = 2;
                    rwName = 'rwY';
                    rwDcm = [1       0             0
                              0   cos(-pi / 2)  -sin(-pi / 2)
                              0   sin(-pi / 2)   cos(-pi / 2)];
                case 'Z'
                    rwIdx = 3;
                    rwName = 'rwZ';
                    rwDcm = eye(3);
            end

            this.reactionwheels(rwIdx).setName(rwName);
            this.reactionwheels(rwIdx).setDcm(rwDcm);
        end
    end
end