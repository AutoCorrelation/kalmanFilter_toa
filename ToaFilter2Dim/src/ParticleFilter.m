classdef ParticleFilter
    properties
        processNoise
        toaNoise
        numParticles
    end

    methods
        function obj = ParticleFilter(Noise, numParticles)
            obj.numParticles = numParticles;
            switch Noise
                case 1
                    obj.processNoise = load('../data/processNoise1.csv');
                    obj.toaNoise = load('../data/toaNoise1.csv');
                case 2
                    obj.processNoise = load('../data/processNoise2.csv');
                    obj.toaNoise = load('../data/toaNoise2.csv');
                case 3
                    obj.processNoise = load('../data/processNoise3.csv');
                    obj.toaNoise = load('../data/toaNoise3.csv');
                case 4
                    obj.processNoise = load('../data/processNoise4.csv');
                    obj.toaNoise = load('../data/toaNoise4.csv');
                case 5
                    obj.processNoise = load('../data/processNoise5.csv');
                    obj.toaNoise = load('../data/toaNoise5.csv');
                otherwise
                    error('Invalid Noise value. Please choose a value between 1 and 5.');
            end
        end

        function y = sampling(obj, x) % CANNOT ENHANCED
            y = zeros(2, obj.numParticles);
            for k = 1:obj.numParticles
                index = ceil(size(obj.toaNoise, 2) * rand);
                % index = k;
                y(:, k) = x + obj.toaNoise(:, index);
            end
        end

        function y = predict(obj, x, B, u)
            y = zeros(size(x));
            for k = 1:obj.numParticles
                index = ceil(size(obj.processNoise, 2) * rand);
                % index = k;
                y(:, k) = x(:, k) + B(:, k) * u + obj.processNoise(:, index);
            end
        end

        function y = predictParam(obj, x, B, u, countStep, gamma)
            y = zeros(size(x));
            for k = 1:obj.numParticles
                index = ceil(size(obj.processNoise, 2) * rand);
                noise = obj.processNoise(:, index);
                % index = k;
                y(:, k) = x(:, k) + B(:, k) * u + noise * exp(-gamma*(countStep-2));
                % y(:, k) = x(:, k) + B(:, k) * u + noise * gamma^(countStep-2);
            end
        end

        function y = update(~, x, w, z, pinvH, R)
            y = zeros(size(w));
            R = R + 1e-6 * eye(size(R));
            for k = 1:length(w)
                y(k) = w(k) * mvnpdf(z, pinvH*x(:, k), R);
                % y(k) = w(k) * mvnpdf(pinvH * z, x(:, k), pinvH * R * pinvH');
            end
            y = y / sum(y);
        end

        function y = updateParam(~, x, w, z, pinvH, R, gamma)
            y = zeros(size(w));
            R = R + 1e-6 * eye(size(R));
            for k = 1:length(w)
                y(k) = w(k) * mvnpdf(z, pinvH*x(:, k), R*gamma);
                % y(k) = w(k) * mvnpdf(pinvH * z, x(:, k), pinvH * R * pinvH');
            end
            y = y / sum(y);
        end

        function y = estimate(~, x, w)
            y = 0;
            for k = 1:length(w)
                y = y + w(k) * x(:, k);
            end
        end

        function y = resample(~, x, w)
            var_accum = 0;
            Npt = length(w);
            for ind = 1:Npt
                var_accum = var_accum + w(ind)^2;
            end
            Ess = 1 / var_accum;
            if Ess < Npt*2 / 3
                wtc = cumsum(w);
                rpt = rand(Npt, 1);
                [~, ind1] = sort([rpt; wtc]);
                ind = find(ind1 <= Npt) - (0:Npt-1)';
                y = x(:, ind);
            else
                y = x;
            end
        end

        function y = resample2(~, x, w)
            var_accum = 0;
            Npt = length(w);
            for ind = 1:Npt
                var_accum = var_accum + w(ind)^2;
            end
            Ess = 1 / var_accum;
            if Ess < Npt*2 / 3
                c = cumsum(w);
                for j = 1:Npt
                    r = rand;
                    b = find(r <= c, 1, 'first');
                    if isempty(b)
                        b = Npt; % 만약 find가 빈 배열을 반환하면 마지막 인덱스를 사용
                    end
                    y(:,j) = x(:,b);
                end
            else
                y = x;
            end
        end
    end
end