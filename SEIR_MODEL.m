I0=0.001*rand(1);
S=1-I0;
E=0;
I=I0;

r1 = 0.1;
r2 = 0.3; 
b1 = 0.005;
b2 = 0.1;
b3 = 0.05;
b4 = 0.2;
dt=1;


for i=1:200
    y(i) = (0.2*E(i)+I(i)) * exp(randn(1)*0.3);
    mul_dist1 = exp(randn(1)*0.5);
    mul_dist2 = exp(randn(1)*0.3);
    mul_dist3 = exp(randn(1)*0.3);

    S(i+1) = S(i)+dt*( (-r1*E(i)*S(i) - r2* I(i) * S(i))*mul_dist1  + b1 * (1-E(i)-I(i)-S(i))*mul_dist2 );
    E(i+1) = E(i)+dt*( (r1*S(i)*E(i) + r2* I(i) * S(i))*mul_dist1- mul_dist3* (b2+b3) * E(i));
    I(i+1) = I(i)+dt*mul_dist3*(b2 * E(i)-b4 *I(i));
    
end


%figure
%plot(S)
%hold on
%plot(I)
%plot(E)
%plot(1-S-I-E)
%plot(y)