classdef MtqArray < handle

    properties(SetAccess = protected, GetAccess = public)
        magnetorquers           % array of magnetorquers (X, Y , and Z unless set up otherwise)
        maxMagneticMoment       % 3d vector with max absolute values of magnetic moment
    end

    methods
        function this = MtqArray(parameters)
            arguments
                parameters.baselineMtq Magnetorquer
                parameters.maxInputVoltage {mustBePositive} = 5;
                parameters.mtqCount {mustBeNumericOrLogical} = 3;
                parameters.doStandardXyzArray {mustBeNumericOrLogical} = true;
            end

            mtqArray(1, parameters.mtqCount) = parameters.baselineMtq;
            this.magnetorquers = mtqArray;

            if parameters.doStandardXyzArray
                for mtqDirection = ['X', 'Y', 'Z']
                    this.updateMagnetorquer(mtqDirection)
                end

                this.updateMaxMagneticMoment(parameters.maxInputVoltage);
            end
        end
    end

    methods(Access = private)

        function updateMagnetorquer(this, mtqDirection)
            switch mtqDirection
                case 'X'
                    mtqIdx = 1;
                    mtqName = 'mtqX';
                    mtqDcm = [cos(pi / 2) 0 sin(pi / 2)
                                  0       1     0
                             -sin(pi / 2) 0 cos(pi / 2)];
                case 'Y'
                    mtqIdx = 2;
                    mtqName = 'mtqY';
                    mtqDcm = [1       0             0
                              0   cos(-pi / 2)  -sin(-pi / 2)
                              0   sin(-pi / 2)   cos(-pi / 2)];
                case 'Z'
                    mtqIdx = 3;
                    mtqName = 'mtqZ';
                    mtqDcm = eye(3);
            end

            this.magnetorquers(mtqIdx).setName(mtqName);
            this.magnetorquers(mtqIdx).setDcm(mtqDcm);
        end

        % this is only valid for standard XYZ configuration, where the
        % directions of output magnetic moment are orthgonal
        function updateMaxMagneticMoment(this, maxVoltage)
            this.maxMagneticMoment = zeros(3, 1);

            for mtqIdx = 1:size(this.magnetorquers, 2)
                this.maxMagneticMoment = this.maxMagneticMoment + this.magnetorquers(mtqIdx).calcMaxMagneticMoment(maxVoltage);
            end
        end
    end
end