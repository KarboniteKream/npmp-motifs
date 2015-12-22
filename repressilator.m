function repressilator()
t = 100;
dt = 0.01;

G = zeros(3,10);

for i = 1:3
    j = i + 1;
    if j == 4
        j = 1;
    end
    
    G(i,1) = 0;
    G(i,2) = 100;
    G(i,3) = -j;
    G(i,4) = 1;
    G(i,5) = 0;
    G(i,6) = 0;
    G(i,7) = 0;
    G(i,8) = 1;
end
setGlobalx(G);

P = zeros(1,3);
P(1) = 100;

figure(1)
[T, P1]=ode15s(@model_complete, [0,t-dt], P);
plot(T,P1)
xlabel('Time')
ylabel('Protein concentration')

end