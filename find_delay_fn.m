function [delay,amp,dc,err]=find_delay_fn(X,Y,f,plot_en)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function [delay,amp,dc,err]=find_delay_fn(X,Y,f,plot_en)
%   Calculate time delay of sensed pulse signal (Sinewave)	
%   Inputs
%   X: time
%   Y: sinwave
%   f: sinwave frequency   
%   plot_en: Enable plot for debugging
%   Outputs
%   dc: DC offset
%   amp: Sinewave amplitude
%   delay: Time delay
%   err: Fitting error

[B0, B1, B2]=sinefit_fn(X,Y,f,plot_en);
delay=(B2+0.5*(B1>0))*-1;
if delay>1/f/2
    delay=mod(delay,1/f);
end
Yhat = B0 + B1*sin(2*pi*f*(X-B2)-pi/2);
err=sqrt(mean((Y-Yhat).^2));
amp=B1;
dc=B0;
