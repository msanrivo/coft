figure (1)
c=load('salida_cowell.txt');
v=load('salida.txt');
t=0:10:86400;
semilogy(t, abs(c(:,1)-v(:,1))/abs(c(:,1)))
%plot(t, c(:,1), t, v(:,1))
grid on
hold on
xlabel('t (s)');
ylabel(' Error relativo en a ');
%ylabel(' a (m) ');
%legend('Cowell', 'VOP')

figure(2)
semilogy(t, abs(c(:,2)-v(:,2))/abs(c(:,2)))
%plot(t, c(:,2), t, v(:,2))
grid on
hold on
xlabel('t (s)');
ylabel(' Error relativo en e ');
%ylabel(' e ');
%legend('Cowell', 'VOP')


figure(3)
semilogy(t, abs(c(:,3)-v(:,3))/abs(c(:,3)))
%plot(t, c(:,3), t, v(:,3))
grid on
hold on
xlabel('t (s)');
ylabel(' Error relativo en i ');
%ylabel(' i (deg) ');
%legend('Cowell', 'VOP')

figure(4)
%plot(t, c(:,4), t, v(:,4))
semilogy(t, abs(c(:,4)-v(:,4))/abs(c(:,4)))
grid on
hold on
xlabel('t (s)');
ylabel(' Error relativo en O ');
%ylabel(' O (deg) ');
%ylabel(' Error relativo en O ');
%legend('Cowell', 'VOP')

figure(5)
%plot(t, c(:,5), t, v(:,5))
semilogy(t, abs(c(:,5)-v(:,5))/abs(c(:,5)))
grid on
hold on
xlabel('t (s)');
ylabel(' Error relativo en w ');
%ylabel(' w (deg) ');
%legend('Cowell', 'VOP')

figure(6)
plot(t, abs(c(:,6)-v(:,6))/abs(c(:,6)))
%plot(t, c(:,6), t, v(:,6))
grid on
hold on
xlabel('t (s)');
%ylabel(' M (deg) ');
ylabel(' Error relativo en M ');
%legend('Cowell', 'VOP')


