function particles = update_location(particles, dt)
    particles(:,1) = particles(:,1) + particles(:,3).*dt;
    particles(:,2) = particles(:,2) + particles(:,4).*dt;
end