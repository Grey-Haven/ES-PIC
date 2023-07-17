function particles = remove_particles(particles, ax, bx, ay, by)
    preserved_particle_idxs = find(particles(:,1) > ax & particles(:,1) < bx ...
                                 & particles(:,2) > ay & particles(:,2) < by);
    particles = particles(preserved_particle_idxs,:);
end