function [xhat, Phat] = ToaKf(x, P, B, u, A, Q, bias, H, R, z)
    % Kalman Filter
    % Prediction
    xhat = A * x + B * u + bias;
    Phat = A * P * A' + Q;
    
    % Estimation
    
    K = Phat * H' * pinv(H * Phat * H' + R);
    xhat = xhat + K * (z - H * xhat);
    Phat = (eye(size(K, 1)) - K * H) * Phat;
end