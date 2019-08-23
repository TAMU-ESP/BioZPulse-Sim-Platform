function [B0h, B1h, B2h]=sinefit_fn(X,Y,f,plot_en)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function [B0h, B1h, B2h]=sinefit_fn(X,Y,f,plot_en)
%   Sinewave paramters extraction by fitting to a regression model 
%   Inputs
%   X: time
%   Y: sinwave
%   f: sinwave frequency   
%   Outputs
%   B0h: DC offset
%   B1h: Sinewave amplitude
%   B2h: Phase shift


B0=mean([max(Y),min(Y)]);
B1 = (max(Y) - min(Y))/2; 
B2 = 0; 
str= strcat('y ~ b0 + b1*sin(2*pi*',num2str(f),'*(x1-b2) - pi/2)');
myFit = NonLinearModel.fit(X,Y,str, [B0, B1, B2]);

B0h= myFit.Coefficients.Estimate(1);
B1h= myFit.Coefficients.Estimate(2);
B2h= myFit.Coefficients.Estimate(3);


if plot_en==1
figure
scatter(X,Y)
grid on;hold on;
Yhat = B0h + B1h*sin(2*pi*f*(X-B2h)-pi/2);
plot(X, Yhat)
end

