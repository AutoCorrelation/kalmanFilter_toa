function [xm, xcov] = UT(Xi, W, noisecov)

    [n, kmax] = size(Xi);
    xm = 0;
    for k = 1:kmax
        xm = xm + W(k)*Xi(:, k);
    end

    xcov = zeros(n, n);
    for k = 1:kmax
        xcov = xcov + W(k)*(Xi(:, k) - xm)*(Xi(:, k) - xm)';
    end
    
    xcov = xcov + noisecov;
end