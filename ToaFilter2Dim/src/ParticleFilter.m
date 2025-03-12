classdef ParticleFilter
    properties
        errorPos
    end

    methods
        function obj = ParticleFilter(Noise)
            switch Noise
                case 1
                    obj.errorPos = load('../data/processNoise1.csv');
                case 2
                    obj.errorPos = load('../data/processNoise2.csv');
                case 3
                    obj.errorPos = load('../data/processNoise3.csv');
                case 4
                    obj.errorPos = load('../data/processNoise4.csv');
                case 5
                    obj.errorPos = load('../data/processNoise5.csv');
                otherwise
                    error('Invalid Noise value. Please choose a value between 1 and 5.');
            end
        end

        function y = predict(x, B, u)
            y = zeros(size(x));
            for k = 1:size(x, 2)
                index = ceil(size(obj.errorPos, 2) * rand);
                y(:, k) = x(:, k) + B * u + obj.errorPos(:, index);
            end
        end

        function y = update(x, w, z, pinvH, R)
            y = zeros(1, size(x, 2));
            for k = 1:size(x, 2)
                y(1, k) = w * mvnpdf(pinvH * z, x, pinvH * R * pinvH');
            end
            y = y / sum(y);
        end

        function y = estimate(x, w)
            y = 0;
            for k = 1:size(x, 2)
                y = y + w(1, k) * x(:, k);
            end
        end
    end
end