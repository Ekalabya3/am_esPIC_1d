% Function for interpolation
function interp = InterpolationF(Interpolation_method, Dx, Ng, position, N_particles, aux_vector)
    switch Interpolation_method
        case 'Nearest Grid Point (NGP)'
            if N_particles == 0
                interp = 0;
            else
                im = floor(position/Dx);
                ip = im + 1;
                Project = [im; ip];
                Project = PeriodicBC(Project, 1, Ng);
                Fim = 1 - abs((position/Dx) - im);
                Fip = 1 - Fim;
                Fraction = [Fim; Fip];
                Fraction(Fraction > 0.5) = 1;
                Fraction(Fraction < 0.5) = 0;
                interp = sparse(aux_vector, Project, Fraction, N_particles, Ng);
            end
        case 'Cloud in Cell (CIC)'
            if N_particles == 0
                interp = 0;
            else
                im = floor(position/Dx);
                ip = im + 1;
                Project = [im; ip];
                Project = PeriodicBC(Project, 1, Ng);
                Fim = 1 - abs((position/Dx) - im);
                Fip = 1 - Fim;
                Fraction = [Fim; Fip];
                interp = sparse(aux_vector, Project, Fraction, N_particles, Ng);
            end
    end
end

