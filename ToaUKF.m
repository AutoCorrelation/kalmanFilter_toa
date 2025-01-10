function [xhat, Phat] = ToaUKF(x, P, B, u, A, Q, bias, H, R, z)
    [Xi, W] = SigmaPoints(x, P, 1);
    fXi = zeros(2, 5);
    for k = 1:5
        fXi(:, k) = A * Xi(:, k) + B * u + bias;
    end

    [xp, Pp] = UT(fXi, W, Q);

    hXi = zeros(6, 5);
    for k = 1:5
        hXi(:, k) = H * fXi(:, k);
    end

    [zp, Pz] = UT(hXi, W, R);

    Pxz = zeros(2, 6);
    for k = 1:5
        Pxz = Pxz + W(k) * (fXi(:, k) - xp) * (hXi(:, k) - zp)';
    end

    K = Pxz * pinv(Pz);

    xhat = xp + K * (z - zp);
    Phat = Pp - K * Pz * K';
end

%-------------------------
function [Xi, W] = SigmaPoints(xm, P, kappa)

    n = numel(xm);
    Xi = zeros(n, 2*n+1);
    W = zeros(2*n+1, 1);

    Xi(:, 1) = xm;
    W(1) = kappa/(n+kappa);

    U = chol((n+kappa)*P);

    for k = 1:n
        Xi(:, k+1) = xm + U(k, :)';
        W(k+1) = 1/(2*(n+kappa));
    end

    for k = 1:n
        Xi(:, n+k+1) = xm - U(k, :)';
        W(n+k+1) = 1/(2*(n+kappa));
    end    
end