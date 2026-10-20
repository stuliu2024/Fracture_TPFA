function c = regular_case21
c.vertices=[
 .5 0 1, .5 0 0, .5 1 0, .5 1 1;
 0 .5 0, 1 .5 0, 1 .5 1, 0 .5 1;
 0 0 .5, 1 0 .5, 1 1 .5, 0 1 .5;
 .75 .5 1, .75 .5 .5, .75 1 .5, .75 1 1;
 .5 .5 .75, 1 .5 .75, 1 1 .75, .5 1 .75;
 .5 .75 .5, 1 .75 .5, 1 .75 1, .5 .75 1;
 .5 .625 .5, .75 .625 .5, .75 .625 .75, .5 .625 .75;
 .625 .5 .75, .625 .5 .5, .625 .75 .5, .625 .75 .75;
 .5 .5 .625, .75 .5 .625, .75 .75 .625, .5 .75 .625];
c.aperture=1e-4;
c.tangential_conductivity=1; % aperture * fracture hydraulic conductivity, m^2/s
c.matrix_high=1;c.matrix_low=.1; % m/s
c.low_boxes=[.5 1 0 .5 0 1; .75 1 .5 .75 .5 1; .625 .75 .5 .625 .5 .75];
c.inlet_interval=[0 .25];c.inlet_flux_density=1; % inward m/s
c.outlet_interval=[.875 1];c.outlet_head=1; % m
c.profile_endpoints=[0 0 0;1 1 1];
end
