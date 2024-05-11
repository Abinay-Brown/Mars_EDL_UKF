function rho = Density(h)
    rho0 = 0.02;
    H = 11.1;               % Mars scale height (km)
    rho = rho0 .* exp(-h./H./1000);  % Atmospheric Density (kg/m^3)
end