function net=grid_gen_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=grid_gen_fn(net)
%   Generate Impedance Grid
%   Outputs: net.grid.x,net.grid.y,net.grid.z

small_edges=0;

Ni=net.Ni;
Nj=net.Nj;
Nk=net.Nk;

% net.param.ri_dc=net.param.r_dc;
% net.param.ci_dc=net.param.c_dc;
net.param.small_edges=small_edges;

% Horizontal Resistors (x)
x=[];
for i=1:Ni+1
    for j=1:Nj
        for k=1:Nk+1
            if net.param.small_edges==1 && (j==1||j==Nj)
                x.r_dc_arr(i,j,k)=net.param.r_dc/2;
            else
                x.re_dc_arr(i,j,k)=([net.re_img(i,j,k)]);
                x.ri_dc_arr(i,j,k)=([net.ri_img(i,j,k)]);
                x.cm_dc_arr(i,j,k)=([net.cm_img(i,j,k)]);
            end
        end
    end
end

% Vert. Resistors (y)
y=[];
for i=1:Ni
    for j=1:Nj+1
        for k=1:Nk+1
            if net.param.small_edges==1 && (i==1||i==Ni)
                y.r_dc_arr(i,j,k)=net.param.r_dc/2;
            else
                y.re_dc_arr(i,j,k)=([net.re_img(i,j,k)]);
                y.ri_dc_arr(i,j,k)=([net.ri_img(i,j,k)]);
                y.cm_dc_arr(i,j,k)=([net.cm_img(i,j,k)]);
            end
        end
    end
end

% Depth Resistors (z)
z=[];
for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk
            if net.param.small_edges==1 && (k==1||k==Nk)
                z.r_dc_arr(i,j,k)=net.param.r_dc/2;
            else
                z.re_dc_arr(i,j,k)=([net.re_img(i,j,k)]);
                z.ri_dc_arr(i,j,k)=([net.ri_img(i,j,k)]);
                z.cm_dc_arr(i,j,k)=([net.cm_img(i,j,k)]);
            end
        end
    end
end


net.grid.x=x;
net.grid.y=y;
net.grid.z=z;



