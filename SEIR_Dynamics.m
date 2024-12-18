function particles = SEIR_Dynamics(S, E, I, w_k1, w_k2, w_k3)
r1 = 0.1;
r2 = 0.3; 
b1 = 0.005;
b2 = 0.1;
b3 = 0.05;
b4 = 0.2;
dt=1;


particles = zeros(3,1);

particles(1) = S - dt*(r1*S*E + r2*S*I)*w_k1 + dt*b1*(1 - S - E - I)*w_k2;
particles(2) = E + dt*(r1*S*E + r2*S*I)*w_k1 - (b2 + b3)*E*w_k2;
particles(3) = I + dt*(b2*E - b4*I)*w_k3;

end

