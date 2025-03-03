function [xhat, pt] = ToaPF(x, B, u, A, Q, bias, H, z)
    % Idea of the Particle Filter
    % Weight dimension: apply weight to x vector
    % Generate particles based on pre simulation.
    persistent firstRun
    Npt = 1e3;
    pt = zeros(2, Npt);
    if isempty(firstRun)
        firstRun = 1;
        % Generate particles based on GAUSSIAN distribution
        pt(1,:) = x(1) + randn(1, Npt);
        pt(2,:) = x(2) + randn(1, Npt);
    end
    
    wt = ones(1, Npt) / Npt;
    pseudoInverseH = inv(H'*H)*H';
    meas = pseudoInverseH*z;
    temp_wt = zeros(2,Npt);
    % Prediction
    for k = 1:Npt
        pt(:, k) = A * pt(:, k) + B * u + bias + sqrtm(Q) * randn(2, 1);
        temp_wt(1,k) = normpdf(meas(1,1), pt(1,k),1); 
        temp_wt(2,k) = normpdf(meas(2,1), pt(2,k),1);
        % wt(1,k) = sqrt(normpdf(meas(1,1), pt(1,k),1) * normpdf(meas(2,1), pt(2,k), 1));
    end
    temp_wt = [temp_wt(1,:)/sum(temp_wt(1,:)); temp_wt(2,:)/sum(temp_wt(2,:))];
    xhat = [pt(1,:)*temp_wt(1,:)'; pt(2,:)*temp_wt(2,:)'];
    % xhat = pt*wt';

    % Resampling
    for i = 1:2
        wtc = cumsum(temp_wt(i,:));
        rpt = rand(Npt,1);
        [~, ind1]= sort([rpt; wtc']);
        ind = find(ind1<=Npt)-(0:Npt-1)';
        pt(i,:) = pt(i,ind);
    end
end

%-------------------------
function normpdf = normpdf(x, mu, sigma)
    normpdf = exp(-0.5*((x-mu)/sigma).^2) / (sigma*sqrt(2*pi));
end