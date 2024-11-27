function value = pdfZDiv(z, options)
    arguments
        z {mustBeNumeric}
        options.Tol (1,1) {mustBeLessThan(options.Tol, 1)} = eps
    end
    
    value = -2 .* (1 + z .* pdfZ(z, Tol=options.Tol));
end