function H = compute_hamiltonian(particles,x,y,phi)

    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    phi_val = gather_field_2D(phi,particles(:,1:2),x,y);
    
    qs = particles(:,5);
    ms = particles(:,6);
    vxs = particles(:,3);
    vys = particles(:,4);

    U = sum(qs./ms .* phi_val);
    T = .5 * sum((vxs.^2 + vys.^2));

    H = T + U;
end