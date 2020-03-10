function [n]=xyz2node_fn(x,y,z,Nx,Ny,Nz)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function [n]=xyz2node_fn(x,y,z,Nx,Ny,Nz)
%   Convert from X,Y,Z coordinates to node index
%   Inputs
%   x: X-axis coordinates
%   y: Y-axis coordinates
%   z: Z-axis coordinates
%   Nx: Number of nodes in X-axis
%   Ny: Number of nodes in Y-axis
%   Nz: Number of nodes in Z-axis
%   Outputs
%   n: Node index


n=(x-1)*(Ny+1)+(y-1)+(z-1)*(Nx+1)*(Ny+1)+1;