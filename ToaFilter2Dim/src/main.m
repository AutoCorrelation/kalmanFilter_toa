clear all;
close all;
clc;
yesorno = input('do preSimulate? Y/N: ','s');
if yesorno == 'Y'
    Env = Env(1e4);
    Env.preSimulate();
end

% load data
load('../data/z.mat');
load('../data/toaPos.mat');
load('../data/R.mat');
%
RSME = RSME();
% parameters
params = struct();
params.numParticles = 1000;
params.numIterations = 100; %size(toaPos, 2);
params.numPoints = size(toaPos, 3);
params.numNoise = size(toaPos, 4);
params.H = [...
    0, -20
    20, -20
    20, 0
    20, 0
    20, 20
    0, 20];
pinvH = pinv(params.H);

%% Particlefilter
pf_data = struct();
pf_data.particles = zeros(2, params.numParticles, params.numPoints, params.numIterations, params.numNoise);
pf_data.vel = zeros(size(pf_data.particles));
pf_data.weights = params.numParticles \ ones(params.numParticles, params.numPoints, params.numIterations, params.numNoise);
pf_data.estimatedPos = zeros(2, params.numPoints, params.numIterations, params.numNoise);
pf_data.RMSE = zeros(params.numNoise, 1);

for countNoise = 1:params.numNoise
% for countNoise = 5
    pf = ParticleFilter(countNoise, params.numParticles);
    
    for countIter = 1:params.numIterations
        for countPoint = 2:params.numPoints
            if countPoint < 3
                pf_data.estimatedPos(:, countPoint-1, countIter, countNoise) = toaPos(:, countIter, countPoint-1, countNoise);
                pf_data.estimatedPos(:, countPoint, countIter, countNoise) = toaPos(:, countIter, countPoint, countNoise);
                pf_data.particles(:, :, countPoint-1, countIter, countNoise) = pf.sampling(toaPos(:, countIter, countPoint-1, countNoise)); % compare rand with countIter
                pf_data.particles(:, :, countPoint, countIter, countNoise) = pf.sampling(toaPos(:, countIter, countPoint, countNoise));
                pf_data.vel(:, :, countPoint, countIter, countNoise) = pf_data.particles(:, :, countPoint, countIter, countNoise) - pf_data.particles(:, :, countPoint-1, countIter, countNoise);
            else
                % pf_data.particles(:, :, countPoint, countIter) = predict(pf, pf_data.particles(:, :, countPoint-1, countIter), pf_data.vel(:, :, countPoint-1, countIter), 1);
                pf_data.particles(:, :, countPoint, countIter, countNoise) = pf.predictParam(pf_data.particles(:, :, countPoint-1, countIter, countNoise), pf_data.vel(:, :, countPoint-1, countIter, countNoise), 1, countPoint, 3);

                pf_data.weights(:, countPoint, countIter, countNoise) = pf.update(pf_data.particles(:, :, countPoint, countIter, countNoise), pf_data.weights(:, countPoint, countIter, countNoise), z(:, countIter, countPoint, countNoise), params.H, R(:, :, countIter, countPoint, countNoise));
                % pf_data.weights(:, countPoint, countIter) = updateParam(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter), z(:, countIter, countPoint, countNoise), params.H, R(:, :, countIter, countPoint, countNoise),0.3);

                pf_data.estimatedPos(:, countPoint, countIter, countNoise) = pf.estimate(pf_data.particles(:, :, countPoint, countIter, countNoise), pf_data.weights(:, countPoint, countIter, countNoise));
                % resampling --------------------------------------------------------
                pf_data.particles(:, :, countPoint, countIter, countNoise) = pf.resample(pf_data.particles(:, :, countPoint, countIter, countNoise), pf_data.weights(:, countPoint, countIter, countNoise));
                % pf_data.particles(:, :, countPoint, countIter) = metropolis_resampling(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter));
                % pf_data.particles(:, :, countPoint, countIter) = multinomial_resampling(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter));
                % pf_data.particles(:, :, countPoint, countIter) = systematic_resampling(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter));
                % pf_data.particles(:, :, countPoint, countIter) = stratified_resampling(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter));
                % pf_data.particles(:, :, countPoint, countIter) = residual_resampling(pf, pf_data.particles(:, :, countPoint, countIter), pf_data.weights(:, countPoint, countIter));
                
                pf_data.vel(:, :, countPoint, countIter, countNoise) = pf_data.estimatedPos(:, countPoint, countIter, countNoise)*ones(1,params.numParticles) - pf_data.particles(:, :, countPoint-1, countIter, countNoise);
            end
        end
    end
end
pf_data = RSME.getRSME(pf_data.estimatedPos);

%% Kalman Filter
kf_data = struct();
kf_data.estimatedPos = zeros(2, params.numPoints, params.numIterations, params.numNoise);
kf_data.errCov = zeros(2, 2, params.numPoints, params.numIterations, params.numNoise);
kf_data.vel = zeros(size(kf_data.estimatedPos));
kf_data.RMSE = zeros(params.numNoise, 1);

for countNoise = 1:params.numNoise
% for countNoise = 5
    kf = KalmanFilter(countNoise, params.H);
    for countIter = 1:params.numIterations
        for countPoint = 2:params.numPoints
            if countPoint < 3
                kf_data.estimatedPos(:, countPoint-1, countIter, countNoise) = toaPos(:, countIter, countPoint-1, countNoise);
                kf_data.estimatedPos(:, countPoint, countIter, countNoise) = toaPos(:, countIter, countPoint, countNoise);
                kf_data.vel(:, countPoint, countIter, countNoise) = kf_data.estimatedPos(:, countPoint, countIter, countNoise) - kf_data.estimatedPos(:, countPoint-1, countIter, countNoise);
            else
                [xhat, Phat] = kf.predict(kf_data.estimatedPos(:, countPoint-1, countIter, countNoise), kf_data.errCov(:, :, countPoint-1, countIter, countNoise), kf_data.vel(:, countPoint-1, countIter, countNoise), 1);
                kf = kf.update(Phat, R(:, :, countIter, countPoint, countNoise));
                [kf_data.estimatedPos(:, countPoint, countIter, countNoise), kf_data.errCov(:,:, countPoint, countIter, countNoise)] = kf.estimate(xhat, Phat, z(:, countIter, countPoint, countNoise));
                kf_data.vel(:, countPoint, countIter, countNoise) = kf_data.estimatedPos(:, countPoint, countIter, countNoise) - kf_data.estimatedPos(:, countPoint-1, countIter, countNoise);
            end
        end
    end
end

kf_data.RMSE = RSME.getRSME(kf_data.estimatedPos);