% Function for charge density calculation
function rho = Charge_density(charge, interp, dx)
    rho = full((charge/dx)*sum(interp))';
end

