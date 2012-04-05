subplot(2,2,1), plot(labels1,hist_n,'-*b');
legend('Distanza nonridotta');
xlabel('Distanza');
ylabel('Frequenza');

subplot(2,2,2), plot(labels2,hist_r,'-*b');
legend('Distanza ridotta');
xlabel('Distanza');
ylabel('Frequeza');

subplot(2,2,3), plot(labels3,hist_t,'-*r');
legend('Distanza di Rohlin topologica');
xlabel('Distanza');
ylabel('Frequenza');

subplot(2,2,4), plot(labels4,hist_t_r,'-*r');
legend('Distanza topologica ridotta');
xlabel('Distanza');
ylabel('Frequeza');

figure


subplot(2,2,1), plot(labels6,hist_f,'-*b');
legend('Distanza tra partizioni con salto 10');
xlabel('Distanza');
ylabel('Frequenza');

subplot(2,2,3), plot(labels8,hist_f_r,'-*b');
legend('Distanza dopo riduzione con fattori dicotomici');
xlabel('Distanza');
ylabel('Frequeza');

subplot(2,2,2), plot(labels7,hist_f_t,'-*r');
legend('Distanza topologica');
xlabel('Distanza');
ylabel('Frequenza');

subplot(2,2,4), plot(labels9,hist_f_r_t,'-*r');
legend('Distanza topologica dopo la riduzione');
xlabel('Distanza');
ylabel('Frequeza');

