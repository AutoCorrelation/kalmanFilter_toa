classdef simulate
    properties
        iteration
        numPoints
        anchorPos
        rangingInfo
        noiseVariance
        errorPos
        d
        z
        H
        pseudoInverseH
        toaPos
        A
        Q
        w_bias
        P0
        R
        liveR
        reducedR
        kf1x
        kf1P
        ukfx
        ukfP
        pfx
    end

    methods
        function obj = simulate(iteration, numPoints)
            clc;

            obj.iteration = iteration;
            obj.numPoints = numPoints;
            obj.noiseVariance = [0.01; 0.1; 1; 10; 100];
            obj.A = eye(2);
        end

        function obj = setsquareAnchorPos(obj, d)
            obj.d = d;
            obj.anchorPos = [0, d; 0, 0; d, 0; d, d];
        end

        function obj = generateRangingInfo(obj)
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        exactPos = [step-1; step-1];
                        for anchor = 1:4
                            obj.rangingInfo(anchor, iter, step, noise) ...
                                = norm(exactPos - obj.anchorPos(anchor, :)') ...
                                + sqrt(obj.noiseVariance(noise)) * randn;
                        end
                    end
                end
            end
        end

        function obj = toa(obj)
            obj.H = [...
                0, -2*obj.d
                2*obj.d, -2*obj.d
                2*obj.d, 0
                2*obj.d, 0
                2*obj.d, 2*obj.d
                0, 2*obj.d];
            obj.pseudoInverseH = (obj.H' * obj.H) \ obj.H';
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        zvec = [obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(2, iter, step, noise)^2 - obj.d^2;
                            obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(3, iter, step, noise)^2;
                            obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + obj.d^2;
                            obj.rangingInfo(2, iter, step, noise)^2 - obj.rangingInfo(3, iter, step, noise)^2 + obj.d^2;
                            obj.rangingInfo(2, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + 2*(obj.d^2);
                            obj.rangingInfo(3, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + obj.d^2];
                        obj.z(:, iter, step, noise) = zvec;
                        obj.toaPos(:, iter, step, noise) = obj.pseudoInverseH * obj.z(:, iter, step, noise);
                    end
                end
            end
        end

        function obj = getQ(obj)
            obj.Q = zeros(2, 2, length(obj.noiseVariance));
            step = 3;
            pos = step * ones(size(obj.toaPos(:, :, step, :)));
            toaVel = (obj.toaPos(:, :, step, :) - obj.toaPos(:, :, step-1, :)) / 0.1;   % pos 2, 1
            obj.errorPos = pos - obj.toaPos(:, :, step, :) - toaVel * 0.1;              % pos

            toaError = pos - obj.toaPos(:, :, step+1, :);
            obj.w_bias = squeeze(mean(obj.errorPos, 2));    % 2x5
            Eerror = squeeze(mean(toaError, 2));
            x = squeeze(obj.errorPos);
            toaError = squeeze(toaError);
            eeT = zeros(2, 2, obj.iteration, length(obj.noiseVariance));
            xxT = zeros(2, 2, obj.iteration, length(obj.noiseVariance));
            for iter = 1:obj.iteration
                for noise = 1:length(obj.noiseVariance)
                    eeT(:,:,iter,noise) = x(:,iter,noise)*x(:,iter,noise)';

                    xxT(:,:,iter,noise) = toaError(:,iter,noise)*toaError(:,iter,noise)';
                end
            end
            EeeT = squeeze(mean(eeT,3));
            ExxT = squeeze(mean(xxT,3));
            for noise = 1:length(obj.noiseVariance)
                obj.Q(:, :, noise) = EeeT(:,:,noise) - obj.w_bias(:, noise) * obj.w_bias(:, noise)';

                obj.P0(:, :, noise) = ExxT(:, :, noise) - Eerror(:, noise) * Eerror(:, noise)';
            end
        end

        function obj = getR(obj)
            obj.R = zeros(6, 6, obj.numPoints, length(obj.noiseVariance));
            for noise = 1:length(obj.noiseVariance)
                for step = 1:obj.numPoints
                    for iter = 1:obj.iteration
                        obj.R(:, :, step, noise) = obj.R(:, :, step, noise) + ...
                            [4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0;...
                            4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                            4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                            -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                            -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                            0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2)
                            ];
                    end
                    obj.R(:, :, step, noise) = obj.R(:, :, step, noise) / obj.iteration;
                end
            end
        end

        function obj = ToaKalmanFilter(obj)
            obj.kf1x = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            obj.kf1P = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        switch step
                            case 1
                                obj.kf1x(:,iter,step,noise) = [0; 0];
                                obj.kf1P(:,:,iter,step,noise) = obj.P0(:, :, noise);
                                velocity = zeros(2, length(obj.noiseVariance));
                            case 2
                                obj.kf1x(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.kf1P(:,:,iter,step,noise) = obj.kf1P(:,:,iter,step-1,noise);
                            case 3
                                obj.kf1x(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.kf1P(:,:,iter,step,noise) = obj.kf1P(:,:,iter,step-1,noise);
                                velocity(:, noise) = (obj.kf1x(:,iter,step,noise) - obj.kf1x(:,iter,step-1,noise)) / 0.1;
                            otherwise
                                obj.liveR(:,:,iter,step,noise) = [4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0;...
                                    4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                                    4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                                    -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                                    -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                                    0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2)
                                    ];
                                % [obj.kf1x(:,iter,step,noise), obj.kf1P(:,:,iter,step,noise)] = ToaKf(obj.kf1x(:,iter,step-1,noise), obj.kf1P(:,:,iter,step-1,noise), 0.1, velocity(:, noise), obj.A, obj.Q(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                [obj.kf1x(:,iter,step,noise), obj.kf1P(:,:,iter,step,noise)] = ToaKf(obj.kf1x(:,iter,step-1,noise), obj.kf1P(:,:,iter,step-1,noise), 0.1, velocity(:, noise), obj.A, obj.Q(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:,:,iter,step,noise), obj.z(:, iter, step, noise));
                                velocity(:, noise) = (obj.kf1x(:,iter,step,noise) - obj.kf1x(:,iter,step-1,noise)) / 0.1;
                        end
                    end
                end
            end
        end

        function obj = optimizeParam1(obj)
            alphamax = 9;
            kfOptix = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            kfOptiP = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            for a = 1:alphamax
                alpha = 0.1 * a;
                for iter = 1:obj.iteration
                    Qbuf = obj.Q;
                    for step = 1:obj.numPoints
                        for noise = 1:length(obj.noiseVariance)
                            switch step
                                case 1
                                    kfOptix(:,iter,step,noise,a) = [0; 0];
                                    kfOptiP(:,:,iter,step,noise,a) = obj.P0(:, :, noise);
                                    velocity = zeros(2, length(obj.noiseVariance));
                                case 2
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                case 3
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                                otherwise
                                    % [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a)] = ToaKf(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                    [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a)] = ToaKf(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:, :,iter, step, noise), obj.z(:, iter, step, noise));
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                            end
                            % Qbuf(:, :, noise) = Qbuf(:, :, noise)*exp(-alpha*(step-3)); % 0.04 + a * 0.01: 8 7 7 5 3 : 0.12 0.11 0.11 0.09 0.07
                            Qbuf(:, :, noise) = Qbuf(:, :, noise)*alpha; % 0.7~0.8
                        end
                    end
                end
            end
            optiRMSE = zeros(length(obj.noiseVariance),alphamax);
            for iter = 1:obj.iteration
                for step = 2:obj.numPoints
                    pos = [step-1; step-1];
                    for noise = 1:length(obj.noiseVariance)
                        for a = 1:alphamax
                            optiRMSE(noise,a) = optiRMSE(noise,a) + norm(pos - kfOptix(:, iter, step, noise,a));
                        end
                    end
                end
            end
            optiRMSE = optiRMSE / (obj.iteration * (obj.numPoints-1));
            [minRMSE, minIdx] = min(optiRMSE, [], 2);
            % figure;
            semilogx(obj.noiseVariance, minRMSE, 'o-','displayName','Q1');
            disp('optimal gamma: ');
            disp(minIdx);
        end

        function obj = optimizeParam2(obj)
            alphamax = 9;
            kfOptix = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            kfOptiP = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            for a = 1:alphamax
                alpha = 0.15 + 0.05*a;
                for iter = 1:obj.iteration
                    Qbuf = obj.Q;
                    for step = 1:obj.numPoints
                        for noise = 1:length(obj.noiseVariance)
                            switch step
                                case 1
                                    kfOptix(:,iter,step,noise,a) = [0; 0];
                                    kfOptiP(:,:,iter,step,noise,a) = obj.P0(:, :, noise);
                                    velocity = zeros(2, length(obj.noiseVariance));
                                case 2
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                case 3
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                                otherwise
                                    Qbuf(:, :, noise) = obj.Q(:,:,noise)*exp(-alpha*(step-3));

                                    % [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a)] = ToaKf(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                    [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a)] = ToaKf(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:, :,iter, step, noise), obj.z(:, iter, step, noise));
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                            end

                        end
                    end
                end
            end
            optiRMSE = zeros(length(obj.noiseVariance),alphamax);
            for iter = 1:obj.iteration
                for step = 2:obj.numPoints
                    pos = [step-1; step-1];
                    for noise = 1:length(obj.noiseVariance)
                        for a = 1:alphamax
                            optiRMSE(noise,a) = optiRMSE(noise,a) + norm(pos - kfOptix(:, iter, step, noise,a));
                        end
                    end
                end
            end
            optiRMSE = optiRMSE / (obj.iteration * (obj.numPoints-1));
            [minRMSE, minIdx] = min(optiRMSE, [], 2);
            % figure;
            semilogx(obj.noiseVariance, minRMSE, 'o-','displayName','Q2_1');
            disp('optimal gamma: ');
            disp(minIdx);
        end

        function obj = optimizeParam3(obj)
            alphamax = 9;
            numVar = length(obj.noiseVariance);
            kfOptix = zeros(2, obj.iteration, obj.numPoints, numVar,alphamax);
            kfOptiP = zeros(2, 2, obj.iteration, obj.numPoints, numVar,alphamax);
            val = zeros(2, 2, numVar);
            for a = 1:alphamax
                alpha = 0.1 * a;
                for iter = 1:obj.iteration
                    Qbuf = obj.Q;
                    for step = 1:obj.numPoints
                        for noise = 1:numVar
                            switch step
                                case 1
                                    kfOptix(:,iter,step,noise,a) = [0; 0];
                                    kfOptiP(:,:,iter,step,noise,a) = obj.P0(:, :, noise);
                                    velocity = zeros(2, numVar);
                                case 2
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                case 3
                                    kfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    kfOptiP(:,:,iter,step,noise,a) = kfOptiP(:,:,iter,step-1,noise,a);
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                                otherwise
                                    % [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a), val(:,:,noise)] = ToaKf_modified(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                    [kfOptix(:,iter,step,noise,a), kfOptiP(:,:,iter,step,noise,a), val(:,:,noise)] = ToaKf_modified(kfOptix(:,iter,step-1,noise,a), kfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:, :, iter, step, noise), obj.z(:, iter, step, noise));
                                    velocity(:, noise) = (kfOptix(:,iter,step,noise,a) - kfOptix(:,iter,step-1,noise,a)) / 0.1;
                                    Qbuf(:, :, noise) = Qbuf(:, :, noise)*alpha + (1-alpha)*val(:,:,noise); % 0.7~0.8
                            end

                        end
                    end
                end
            end
            optiRMSE = zeros(numVar,alphamax);
            for iter = 1:obj.iteration
                for step = 2:obj.numPoints
                    pos = [step-1; step-1];
                    for noise = 1:numVar
                        for a = 1:alphamax
                            optiRMSE(noise,a) = optiRMSE(noise,a) + norm(pos - kfOptix(:, iter, step, noise,a));
                        end
                    end
                end
            end
            optiRMSE = optiRMSE / (obj.iteration * (obj.numPoints-1));
            [minRMSE, minIdx] = min(optiRMSE, [], 2);
            % figure;
            semilogx(obj.noiseVariance, minRMSE, 'o-','displayName','Q3');
            disp('optimal gamma: ');
            disp(minIdx);
        end

        function obj = ToaUKF(obj)
            obj.ukfx = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            obj.ukfP = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        switch step
                            case 1
                                obj.ukfx(:,iter,step,noise) = [0; 0];
                                obj.ukfP(:,:,iter,step,noise) = obj.P0(:, :, noise);
                                velocity = zeros(2, length(obj.noiseVariance));
                            case 2
                                obj.ukfx(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.ukfP(:,:,iter,step,noise) = obj.ukfP(:,:,iter,step-1,noise);
                            case 3
                                obj.ukfx(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.ukfP(:,:,iter,step,noise) = obj.ukfP(:,:,iter,step-1,noise);
                                velocity(:, noise) = (obj.ukfx(:,iter,step,noise) - obj.ukfx(:,iter,step-1,noise)) / 0.1;
                            otherwise
                                % [obj.ukfx(:,iter,step,noise), obj.ukfP(:,:,iter,step,noise)] = ToaUKF(obj.ukfx(:,iter,step-1,noise), obj.ukfP(:,:,iter,step-1,noise), 0.1, velocity(:, noise), obj.A, obj.Q(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                [obj.ukfx(:,iter,step,noise), obj.ukfP(:,:,iter,step,noise)] = ToaUKF(obj.ukfx(:,iter,step-1,noise), obj.ukfP(:,:,iter,step-1,noise), 0.1, velocity(:, noise), obj.A, obj.Q(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:, :, iter, step, noise), obj.z(:, iter, step, noise));
                                velocity(:, noise) = (obj.ukfx(:,iter,step,noise) - obj.ukfx(:,iter,step-1,noise)) / 0.1;
                        end
                    end
                end
            end
        end

        function obj = optimizeParam4(obj)
            alphamax = 9;
            ukfOptix = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            ukfOptiP = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphamax);
            for a = 1:alphamax
                alpha = 0.1 * a; % 0.45 45 45 35 25
                for iter = 1:obj.iteration
                    Qbuf = obj.Q;
                    for step = 1:obj.numPoints
                        for noise = 1:length(obj.noiseVariance)
                            switch step
                                case 1
                                    ukfOptix(:,iter,step,noise,a) = [0; 0];
                                    ukfOptiP(:,:,iter,step,noise,a) = obj.P0(:, :, noise);
                                    velocity = zeros(2, length(obj.noiseVariance));
                                case 2
                                    ukfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    ukfOptiP(:,:,iter,step,noise,a) = ukfOptiP(:,:,iter,step-1,noise,a);
                                case 3
                                    ukfOptix(:,iter,step,noise,a) = obj.toaPos(:, iter, step, noise);
                                    ukfOptiP(:,:,iter,step,noise,a) = ukfOptiP(:,:,iter,step-1,noise,a);
                                    velocity(:, noise) = (ukfOptix(:,iter,step,noise,a) - ukfOptix(:,iter,step-1,noise,a)) / 0.1;
                                otherwise
                                    Qbuf(:, :, noise) = obj.Q(:,:,noise)*exp(-alpha*(step-3));

                                    % [ukfOptix(:,iter,step,noise,a), ukfOptiP(:,:,iter,step,noise,a)] = ToaUKF(ukfOptix(:,iter,step-1,noise,a), ukfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                    [ukfOptix(:,iter,step,noise,a), ukfOptiP(:,:,iter,step,noise,a)] = ToaUKF(ukfOptix(:,iter,step-1,noise,a), ukfOptiP(:,:,iter,step-1,noise,a), 0.1 , velocity(:, noise), obj.A, Qbuf(:, :, noise), obj.w_bias(:, noise), obj.H, obj.liveR(:, :, iter, step, noise), obj.z(:, iter, step, noise));
                                    velocity(:, noise) = (ukfOptix(:,iter,step,noise,a) - ukfOptix(:,iter,step-1,noise,a)) / 0.1;
                            end

                        end
                    end
                end
            end
            optiRMSE = zeros(length(obj.noiseVariance),alphamax);
            for iter = 1:obj.iteration
                for step = 2:obj.numPoints
                    pos = [step-1; step-1];
                    for noise = 1:length(obj.noiseVariance)
                        for a = 1:alphamax
                            optiRMSE(noise,a) = optiRMSE(noise,a) + norm(pos - ukfOptix(:, iter, step, noise,a));
                        end
                    end
                end
            end
            optiRMSE = optiRMSE / (obj.iteration * (obj.numPoints-1));
            [minRMSE, minIdx] = min(optiRMSE, [], 2);
            % figure;
            semilogx(obj.noiseVariance, minRMSE, 'o-','displayName','Q4');
            disp('optimal gamma: ');
            disp(minIdx);
        end

        function obj = ToaParticleFilter(obj)
            % https://github.com/vatsl/ParticleFilter/blob/master/src/particle_filter.cpp
            totalParticles = 1000;
            obj.pfx = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            for countIter = 1:obj.iteration
                particleTensor = zeros(2, totalParticles, obj.numPoints, length(obj.noiseVariance));
                velocity = zeros(2, totalParticles, obj.numPoints,length(obj.noiseVariance));
                for countStep = 1:obj.numPoints
                    for countNoise = 1:length(obj.noiseVariance)
                        switch countStep
                            case 1
                                obj.pfx(:,countIter, countStep, countNoise) = [0; 0];
                            case 2
                                obj.pfx(:, countIter, countStep, countNoise) = obj.toaPos(:, countIter, countStep, countNoise);
                                particleTensor(:,:,countStep,countNoise) = ...
                                    obj.pfx(:, countIter, countStep, countNoise) + mvnrnd(zeros(2,1), obj.Q(:, :, 5), totalParticles)';
                            case 3
                                obj.pfx(:, countIter, countStep, countNoise) = obj.toaPos(:, countIter, countStep, countNoise);

                                particleTensor(:,:,countStep,countNoise) = ...
                                    obj.pfx(:, countIter, countStep, countNoise) + mvnrnd(zeros(2,1), obj.Q(:, :, 5), totalParticles)';
                                velocity(:, :, countStep,countNoise) = (particleTensor(:,:,countStep,countNoise) - particleTensor(:,:,countStep-1,countNoise)) / 0.1;
                                weight = ones(2, totalParticles) / totalParticles;
                            otherwise
                                obj.reducedR(:,:,countIter, countStep, countNoise) = ...
                                    7 * obj.pseudoInverseH*obj.liveR(:,:,countIter, countStep, countNoise)*obj.pseudoInverseH';

                                [obj.pfx(:, countIter, countStep, countNoise), particleTensor(:,:,countStep,countNoise)]= ...
                                    ToaPF(particleTensor(:,:,countStep-1,countNoise), weight, velocity(:,:,countStep, countNoise), 0.1, obj.A, obj.Q(:,:,countNoise), obj.w_bias(:, countNoise), obj.pseudoInverseH, obj.reducedR(:,:,countIter, countStep, countNoise), obj.z(:, countIter, countStep, countNoise), obj.errorPos(:,:,1,countNoise));
                                velocity(:, :, countStep,countNoise) = (particleTensor(:,:,countStep,countNoise) - particleTensor(:,:,countStep-1,countNoise)) / 0.1;
                        end
                    end
                end
            end
        end


        function obj = optimizeToaParticleFilter(obj)
            totalParticles = 100;
            alphaMax = 9;
            obj.pfx = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance),alphaMax);

            for countAlpha = 1:alphaMax
                for countIter = 1:obj.iteration
                    particleTensor = zeros(2, totalParticles, obj.numPoints, length(obj.noiseVariance),alphaMax);
                    velocity = zeros(2, totalParticles, obj.numPoints,length(obj.noiseVariance),alphaMax);
                    for countStep = 1:obj.numPoints
                        for countNoise = 1:length(obj.noiseVariance)
                            switch countStep
                                case 1
                                    obj.pfx(:,countIter, countStep, countNoise,countAlpha) = [0; 0];
                                case 2
                                    obj.pfx(:,countIter, countStep, countNoise,countAlpha) = obj.toaPos(:, countIter, countStep, countNoise);
                                    particleTensor(:,:,countStep,countNoise) = ...
                                        obj.pfx(:,countIter, countStep, countNoise,countAlpha) + mvnrnd(zeros(2,1), obj.Q(:, :, countNoise), totalParticles)';
                                case 3
                                    obj.pfx(:,countIter, countStep, countNoise,countAlpha) = obj.toaPos(:, countIter, countStep, countNoise);

                                    particleTensor(:,:,countStep,countNoise) = ...
                                        obj.pfx(:,countIter, countStep, countNoise,countAlpha) + mvnrnd(zeros(2,1), obj.Q(:, :, countNoise), totalParticles)';
                                    velocity(:, :, countStep,countNoise) = (particleTensor(:,:,countStep,countNoise) - particleTensor(:,:,countStep-1,countNoise)) / 0.1;
                                    weight = ones(2, totalParticles) / totalParticles;
                                otherwise
                                    obj.reducedR(:,:,countIter, countStep, countNoise) = ...
                                        countAlpha*obj.pseudoInverseH*obj.liveR(:,:,countIter, countStep, countNoise)*obj.pseudoInverseH';

                                    [obj.pfx(:,countIter, countStep, countNoise,countAlpha), particleTensor(:,:,countStep,countNoise)]= ...
                                        ToaPF(particleTensor(:,:,countStep-1,countNoise), weight, velocity(:,:,countStep, countNoise), 0.1, obj.A, obj.Q(:,:,countNoise), obj.w_bias(:, countNoise), obj.pseudoInverseH, obj.reducedR(:,:,countIter, countStep, countNoise), obj.z(:, countIter, countStep, countNoise), obj.errorPos(:,:,1,countNoise));
                                    velocity(:, :, countStep,countNoise) = (particleTensor(:,:,countStep,countNoise) - particleTensor(:,:,countStep-1,countNoise)) / 0.1;
                            end
                        end
                    end
                end
            end

            optiRMSE = zeros(length(obj.noiseVariance),alphaMax);
            for countIter = 1:obj.iteration
                for countStep = 2:obj.numPoints
                    pos = [countStep-1; countStep-1];
                    for countNoise = 1:length(obj.noiseVariance)
                        for countAlpha = 1:alphaMax
                            optiRMSE(countNoise,countAlpha) = optiRMSE(countNoise,countAlpha) + norm(pos - obj.pfx(:, countIter, countStep, countNoise,countAlpha));
                        end
                    end
                end
            end
            optiRMSE = optiRMSE / (obj.iteration * (obj.numPoints-1));
            [minRMSE, minIdx] = min(optiRMSE, [], 2);
            % figure;
            semilogx(obj.noiseVariance, minRMSE, 'o-','displayName','OptiPF');
            disp('optimal gamma: ');
            disp(minIdx);
        end


        function obj = resultPlot(obj)
            toaRMSE = zeros(length(obj.noiseVariance),1);
            kf1RMSE = zeros(length(obj.noiseVariance),1);
            ukfRMSE = zeros(length(obj.noiseVariance),1);
            pfRMSE = zeros(length(obj.noiseVariance),1);
            for countIter = 1:obj.iteration
                for countStep = 2:obj.numPoints
                    pos = [countStep-1; countStep-1];
                    for countNoise = 1:length(obj.noiseVariance)
                        toaRMSE(countNoise,1) = toaRMSE(countNoise,1) + norm(pos - obj.toaPos(:, countIter, countStep, countNoise));
                        kf1RMSE(countNoise,1) = kf1RMSE(countNoise,1) + norm(pos - obj.kf1x(:, countIter, countStep, countNoise));
                        % ukfRMSE(noise,1) = ukfRMSE(noise,1) + norm(pos - obj.ukfx(:, countIter, countStep, countNoise));
                        pfRMSE(countNoise,1) = pfRMSE(countNoise,1) + norm(pos - obj.pfx(:, countIter, countStep, countNoise));

                    end
                end
            end

            % for countStep = 2:obj.numPoints
            %     pos = [countStep-1; countStep-1];
            %     for countNoise = 1:length(obj.noiseVariance)
            %         pfRMSE(countNoise,1) = pfRMSE(countNoise,1) + norm(pos - obj.pfx(:, 1, countStep, countNoise));
            %     end
            % end
            toaRMSE = toaRMSE / (obj.iteration * (obj.numPoints-1));
            kf1RMSE = kf1RMSE / (obj.iteration * (obj.numPoints-1));
            % ukfRMSE = ukfRMSE / (obj.iteration * (obj.numPoints-1));
            pfRMSE = pfRMSE / (obj.iteration * (obj.numPoints-1));

            figure;
            semilogx(obj.noiseVariance, toaRMSE, 'o-','DisplayName','TOA');
            hold on;
            semilogx(obj.noiseVariance, kf1RMSE, 'o-','DisplayName','KF1');
            % semilogx(obj.noiseVariance, ukfRMSE, 'o-','DisplayName','UKF');
            semilogx(obj.noiseVariance, pfRMSE, 'o-','DisplayName','PF');
            xlabel('Noise Variance');
            ylabel('RMSE');
            legend;
        end
    end
end



