classdef KalmanFilter
    %KALMANFILTER
    %   Kalman Filtering for 2D TOA
    %   This class implements the Kalman filter algorithm for 2D Time of Arrival (TOA) localization.
    %   구조체를 이용해서 다른 변수들도 전달할 수 있도록 하고, 특히 Q를 구조체로 전달할 수 있도록 한다.
    %   properties에는 아무것도 넣지 않음.
    
    properties
        pNoiseCov   % Process noise covariance
    end
    
    methods
        function obj = KalmanFilter(Noise)
            %KALMANFILTER 이 클래스의 인스턴스 생성
            %   자세한 설명 위치
            switch Noise
                case 1
                    obj.pNoiseCov = load('../data/Q1.csv');
                case 2
                    obj.pNoiseCov = load('../data/Q2.csv');
                case 3
                    obj.pNoiseCov = load('../data/Q3.csv');
                case 4
                    obj.pNoiseCov = load('../data/Q4.csv');
                case 5
                    obj.pNoiseCov = load('../data/Q5.csv');
                otherwise
                    error('Invalid Noise value. Please choose a value between 1 and 5.');
            end
            
        end
        
        function y = predict(obj,stateStruct, B, u, bias)
            x = stateStruct.state;
            P = stateStruct.errCov;

            y.state = x + B * u + bias;
            y.errCov = P + obj.pNoiseCov;
        end
        

        function y = update(~,stateStruct, R)
            P = stateStruct.errCov;
            R = R + 1e-6 * eye(size(R));
            K = P * H' / (H * P * H' + R);
            y = K;
        end

        function y = estimate(~,stateStruct, K, z, H)
            x = stateStruct.state;
            P = stateStruct.errCov;
            y.state = x + K * (z - H * x);
            y.errCov = (eye(size(P)) - K * H) * P;
        end
    end
end

