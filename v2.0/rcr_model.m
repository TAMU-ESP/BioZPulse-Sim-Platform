function yout = rcr_model(b,xdata)
%b(1)=RE, b(2)=CM, b(3)=RI
yout = zeros(length(xdata),3); % allocate yout

% yout(:,1) = real(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));
% yout(:,2) = imag(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));
% yout(:,3) = abs(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));

yout(:,1) = real(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));
yout(:,2) = imag(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));
yout(:,3) = abs(b(1)*(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))./(b(1)+(b(3)+1./(1i*2*pi*xdata*b(2)*1e-9))));
