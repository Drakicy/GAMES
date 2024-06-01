function value = pdfZ(z, options)
    arguments
        z {mustBeNumeric}
        options.Tol (1,1) {mustBeLessThan(options.Tol, 1)} = eps
    end

    value = zeros(size(z));
    
    ind = (real(z) == 0);
    value(ind) = exp(-imag(z(ind)).^2) .* erfc(imag(z(ind)));
    
    ind = ~ind;
    indNeg = ind & (imag(z) < 0);
    z(indNeg) = conj(z(indNeg));

    N = ceil(log(options.Tol / 2^(3/2)) / log(sqrt(2) - 1)) - 1;
    L = sqrt(N / sqrt(2));
    
    M = 2 * N;
    k = (-M+1:1:M-1)';

    t = L * tan(k * pi / (2 * M));
    f = (L^2 + t.^2) .* exp(-t.^2);
    f = [0; f];
    
    a = real(fft(fftshift(f))) / (2 * M);
    a = flipud(a(2:N+1));
    Z = (L + 1i .* z(ind)) ./ (L - 1i .* z(ind));
    p = polyval(a, Z);
    value(ind) = 2 * p ./ (L - 1i .* z(ind)).^2 + 1 / sqrt(pi) ./ (L - 1i .* z(ind));
    
    value(indNeg) = conj(2 * exp(-z(indNeg).^2) - value(indNeg));
    
    value = sqrt(pi) .* 1i .* value;
end