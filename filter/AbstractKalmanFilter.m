classdef AbstractKalmanFilter < handle

    properties (SetAccess = protected, GetAccess = public)
        odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

        sat % satellite object
        orb % orbit object

        P % error covariance matrix
        R % sensor noise covariance matrix
        Q % process noise covariance matrix

        traceP % trace of error covariance matrix
    end

    methods (Abstract)

        estimate(this)

    end

    methods

        function val = get.traceP(this)
            val = trace(this.P);
        end

    end

    methods (Abstract, Access = protected)

        prediction(this)
        correction(this)

        initProcessCovariance(this)
        initErrorCovariance(this)
        initMeasurementsCovariance(this)

    end

end