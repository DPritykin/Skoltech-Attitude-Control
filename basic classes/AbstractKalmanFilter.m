classdef AbstractKalmanFilter < handle

    properties (SetAccess = protected, GetAccess = public)
        odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

        sat % satellite object
        orb % orbit object

        P % error covariance matrix
        R % sensor noise covariance matrix
        Q % process noise covariance matrix
    end

    methods (Abstract)

        estimate(this)

    end

    methods (Abstract, Access = protected)

        prediction(this)
        correction(this)

        initProcessCovariance(this)
        initErrorCovariance(this)
        initMeasurementsCovariance(this)
    end

end