clear all; close all; clc;

Qsim = linspace(1e5,5e5,4);
Qcut = 8e5;
Q = min(Qsim,Qcut)/Qcut;
V = linspace(0.4,1.2,100);
figure;
hold on
cmap = hsv(length(Q))./1.5;
for i = 1:length(Q)
F = Q(i)./V;
Fsq = (Q(i)./V).^2;
plot(V,F,'color',cmap(i,:),'linewidth',2);
plot(V,Fsq,'--','color',cmap(i,:),'linewidth',2);
end

grid on
box on
xlabel('V_{mode}/(\lambda/n)^3')
ylabel('Fitness value')

axis tight
xlim([0.4 1.2])