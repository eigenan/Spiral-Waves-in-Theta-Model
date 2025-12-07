load("meta_data.mat");

% plot wavenumber-mu graph
figure(1)
plot(mu_list,wavenumber,'b.-','LineWidth',2)
title("$\mu$ vs the wavenumber","interpreter","latex","FontSize",24)
xlabel("$\mu$","interpreter","latex","FontSize",24)
ylabel("wavenumber","interpreter","latex","FontSize",24)

% plot omega-mu graph
figure(2)
plot(mu_list,omega_list,'r.-','LineWidth',2)
title("$\mu$ vs $\omega$","interpreter","latex","FontSize",24)
xlabel("$\mu$","interpreter","latex","FontSize",24)
ylabel("$\omega$","interpreter","latex","FontSize",24)