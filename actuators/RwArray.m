classdef RwArray < handle

    properties(SetAccess = protected, GetAccess = public)
        reactionwheels           % array of reaction wheels(X, Y , and Z unless set up otherwise)
    end

    methods
        function this = RwArray(parameters)
            arguments
                parameters.baselineRw ReactionWheel
                parameters.rwCount {mustBeNumericOrLogical} = 3;
                parameters.doStandardXyzArray {mustBeNumericOrLogical} = true;
            end

            rwArray(1, parameters.rwCount) = ReactionWheel();
            for rwIdx = 1:parameters.rwCount
                rwArray(rwIdx).setMaxTorque(parameters.baselineRw.maxTorque);
                rwArray(rwIdx).setMaxAngMomentum(parameters.baselineRw.maxAngMomentum);
            end

            this.reactionwheels = rwArray;

            if parameters.doStandardXyzArray
                for rwAxisBf = ['X', 'Y', 'Z']
                    this.updateReactionWheel(rwAxisBf)
                end
            end
        end

        function actuatedTorque = actuateCommand(this, commandTorque, duration, rwAngMomentum)
            actuatedTorque = zeros(3, 1);

            for rwIdx = 1:size(this.reactionwheels, 2)
                actuatedTorque = actuatedTorque + ...
                                 this.reactionwheels(rwIdx).actuateControlTorque(commandTorque(rwIdx), ...
                                                                                 duration, ...
                                                                                 rwAngMomentum(rwIdx));
            end
        end
    end

    methods(Access = private)

        function updateReactionWheel(this, rwAxisBf)
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