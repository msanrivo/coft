function MapView(t,X,tt,XX)

[lambda,phi,r] = cart2sph(X(:,1),X(:,2),X(:,3));
load coast
geoshow(lat,long,'Color','black','LineWidth',3);
hold on
plot(lambda*180/pi,phi*180/pi,'.-')
hold on
[lambda2,phi2,r2] = cart2sph(XX(:,1),XX(:,2),XX(:,3));
plot(lambda2*180/pi,phi2*180/pi,'ro','MarkerSize',5);
% for i = 1 : length(t)
%     
%     
%     hold on
%     plot(lambda(1:i),phi(1:i),'.--');
%     
%     title(sprintf('Time t = %.5E [secs]',t(i)));
%     drawnow
%     
% end


end

