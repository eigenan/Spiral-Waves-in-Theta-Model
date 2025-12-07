load("dispersion_data.mat");
figure(3)
plot(mu_list(1:end),omega_list(1:end),'b.-')
hold on
plot(mu_list(1:end),-mu_list(1:end).^2./2+1,'ro')
xlabel('$\mu$','Interpreter','Latex',"FontSize",18)
ylabel('$\omega$','Interpreter','Latex',"FontSize",18)