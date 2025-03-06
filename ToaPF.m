function [xhat, particleTensor] = ToaPF(x, weight, B, u, A, Q, bias, pseudoInverseH, R, z)
    % Idea of the Particle Filter
    % Weight dimension: apply weight to x vector
    % Generate particles based on pre simulation.
    % Generate particles based on GAUSSIAN distribution
    % 초기 난수는 가우시안 분포를 따르는 난수를 사용한다.

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
    % SIR
    % c = cumsum(tempWeight);
    % for j = 1:Npt
    %     r = rand;
    %     b = find(r <= c, 1, 'first');
    %     if isempty(b)
    %         b = Npt; % 만약 find가 빈 배열을 반환하면 마지막 인덱스를 사용
    %     end
    %     particleTensor(:,j) = particleTensor(:,b);
    % end

    wtc = cumsum(tempWeight);
    rpt = rand(Npt,1);
    [~, ind1]= sort([rpt; wtc]);
    ind = find(ind1<=Npt)-(0:Npt-1)';
    particleTensor = particleTensor(:,ind);
end

%-------------------------
% function normpdf = normpdf(x, mu, sigma)
%     normpdf = exp(-0.5*((x-mu)/sigma).^2) / (sigma*sqrt(2*pi));
% end