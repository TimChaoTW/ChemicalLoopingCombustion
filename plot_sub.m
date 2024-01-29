function plot_sub(M_title,y_name,save_name)
title (M_title)
xlabel ('Time sec')
ylabel (y_name)
legend('Z=0.4m','Z=0.8m','Z=1.2m','Z=1.6m','Z=2m','Z=2.2m')
saveas(gcf,save_name,'png')
% plot_end = 0;