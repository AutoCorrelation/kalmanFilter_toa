function [xhat, particleTensor] = ToaPF(x, weight, B, u, A, Q, bias, pseudoInverseH, R, z)
    % Idea of the Particle Filter
    % Weight dimension: apply weight to x vector
    % Generate particles based on pre simulation.
    % Generate particles based on GAUSSIAN distribution
    % 초기 난수는 가우시안 분포를 따르는 난수를 사용한다.
    % P0 = 0.005309	-0.000024
    %       -0.000024	0.005173
    Npt = size(weight, 2);
    particleTensor = zeros(2, Npt);
    xhat = zeros(2, 1);
    tempWeight = zeros(size(weight,2),1);
    reducedZ = pseudoInverseH*z;
    % Prediction
    % 예측 오차랜덤 값을 Q를 기반으로 생성하되 스텝에 따라서 감소하도록 한다. (fading memory Q)
    for k = 1:Npt
        particleTensor(:, k) = A * x(:,k) + B * u + mvnrnd(bias, Q, 1)';  % different random value for each particle(w_k)
        % xhat = xhat + weight(:,k) .* particleTensor(:,k);
        tempWeight(k,1) =  mvnpdf(particleTensor(:,k),reducedZ,R);
    end
    tempWeight = tempWeight / sum(tempWeight);
    for j = 1:Npt
        xhat = xhat + tempWeight(j,1) * particleTensor(:,j);
    end
    % Resampling
    % 중요한 파티클을 더 많이 살리기 위해 중요도 샘플링을 사용한다.
    c = cumsum(tempWeight);
    for j = 1:Npt
        r = rand;
        b = find(r <= c, 1, 'first');
        particleTensor(:,j) = particleTensor(:,b);
    end
end

%-------------------------
% function normpdf = normpdf(x, mu, sigma)
%     normpdf = exp(-0.5*((x-mu)/sigma).^2) / (sigma*sqrt(2*pi));
% end