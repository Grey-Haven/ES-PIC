function particles = update_velocity(particles, x, y, phi, dx, dy, dt)
    phi_x = zeros(size(phi));
    phi_y = zeros(size(phi));
    phi_x(:,1) = (-3*phi(:,1) + 4*phi(:,2) - phi(:,3))/(2*dx);
    phi_x(:,2:end-1) = (phi(:,3:end) - phi(:,1:end-2))/(2*dx);
    phi_x(:,end) = (3*phi(:,end) - 4*phi(:,end-1) + phi(:,end-2))/(2*dx);
    
    
    phi_y(1,:) = (-3*phi(1,:) + 4*phi(2,:) - phi(3,:))/(2*dy);
    phi_y(2:end-1,:) = (phi(3:end,:) - phi(1:end-2,:))/(2*dy);
    phi_y(end,:) = (3*phi(end,:) - 4*phi(end-1,:) + phi(end-2,:))/(2*dy);
    
    Ex = -phi_x;
    Ey = -phi_y;
    
    Ex_vals = gather_field_2D(Ex,particles(:,1:2),x,y);
    Ey_vals = gather_field_2D(Ey,particles(:,1:2),x,y);
    
    qs = particles(:,5);
    ms = particles(:,6);
    
    particles(:,3) = particles(:,3) + qs./ms .* Ex_vals * dt;
    particles(:,4) = particles(:,4) + qs./ms .* Ey_vals * dt;
end