clear all;
close all;
clc;
yesorno = input('do preSimulate? Y/N: ','s');
if yesorno == 'Y'
    Env = Env(1e4);
    Env.preSimulate();
end
%%
% load data
load('../data/z.mat');
load('../data/toaPos.mat');
load('../data/R.mat');
%
% parameters
params = struct();
params.numParticles = 2000;
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

% Particlefilter
RSME = RSME();
pf_RMSE = zeros(params.numNoise, 1);
for countNoise = 1:params.numNoise
% for countNoise = 5
    pf = ParticleFilter(countNoise, params.numParticles);
    pf_particles = zeros(2, params.numParticles, params.numPoints, params.numIterations);
    pf_vel = zeros(size(pf_particles));
    pf_weights = params.numParticles\ones(params.numParticles, params.numPoints, params.numIterations);
    pf_estimatedPos = zeros(2, params.numPoints, params.numIterations);
    pf_errorPos = zeros(size(pf_estimatedPos));
    
    for countIter = 1:params.numIterations
        for countPoint = 2:params.numPoints
            if countPoint < 3
                pf_estimatedPos(:, countPoint-1, countIter) = toaPos(:, countIter, countPoint-1, countNoise);
                pf_estimatedPos(:, countPoint, countIter) = toaPos(:, countIter, countPoint, countNoise);
                pf_particles(:, :, countPoint-1, countIter) = sampling(pf, toaPos(:, countIter, countPoint-1, countNoise)); % compare rand with countIter
                pf_particles(:, :, countPoint, countIter) = sampling(pf, toaPos(:, countIter, countPoint, countNoise));
                pf_vel(:, :, countPoint, countIter) = pf_particles(:, :, countPoint, countIter) - pf_particles(:, :, countPoint-1, countIter);
            else
                % pf_particles(:, :, countPoint, countIter) = predict(pf, pf_particles(:, :, countPoint-1, countIter), pf_vel(:, :, countPoint-1, countIter), 1);
                pf_particles(:, :, countPoint, countIter) = predictParam(pf, pf_particles(:, :, countPoint-1, countIter), pf_vel(:, :, countPoint-1, countIter), 1, countPoint, 0.5);
                pf_weights(:, countPoint, countIter) = update(pf, pf_particles(:, :, countPoint, countIter), pf_weights(:, countPoint, countIter), z(:, countIter, countPoint, countNoise), params.H, R(:, :, countIter, countPoint, countNoise));
                % pf_weights(:, countPoint, countIter) = updateParam(pf, pf_particles(:, :, countPoint, countIter), pf_weights(:, countPoint, countIter), z(:, countIter, countPoint, countNoise), params.H, R(:, :, countIter, countPoint, countNoise),0.3);
                pf_estimatedPos(:, countPoint, countIter) = estimate(pf, pf_particles(:, :, countPoint, countIter), pf_weights(:, countPoint, countIter));
                pf_particles(:, :, countPoint, countIter) = resample(pf, pf_particles(:, :, countPoint, countIter), pf_weights(:, countPoint, countIter));
                % pf_particles(:, :, countPoint, countIter) = resample2(pf, pf_particles(:, :, countPoint, countIter), pf_weights(:, countPoint, countIter));
                pf_vel(:, :, countPoint, countIter) = pf_particles(:, :, countPoint, countIter) - pf_particles(:, :, countPoint-1, countIter);
            end
        end
    end

    % error
    pf_RMSE(countNoise) = getRSME(RSME, pf_estimatedPos);
end

