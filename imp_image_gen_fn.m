function r_dc_arr=imp_image_gen_fn(r_in,net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function r_dc_arr=imp_image_gen_fn(r_in,net)
%   Adjust bio-impoedance values to MATLAB

Ni=net.Ni;
Nj=net.Nj;
Nk=net.Nk;

if length(r_in)==1
    for i=1:Ni+1
        for j=1:Nj+1
            for k=1:Nk+1
                r_img(i,j,k)=r_in;
            end
        end
    end    
else
    r_img=r_in;
end

for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk+1
            r_dc_arr(i,j,k).x=NaN;
            r_dc_arr(i,j,k).y=NaN;
            r_dc_arr(i,j,k).z=NaN;
        end
    end
end

% Horizontal Resistors (x)
for i=1:Ni+1
    for j=1:Nj
        for k=1:Nk+1
            if net.param.small_edges==1 && (j==1||j==Nj)
                r_dc_arr(i,j,k).x=net.param.r_dc/2;
            else
                r_dc_arr(i,j,k).x=r_img(i,j,k);
            end
        end
    end
end

% Vert. Resistors (y)
for i=1:Ni
    for j=1:Nj+1
        for k=1:Nk+1
            if net.param.small_edges==1 && (i==1||i==Ni)
                r_dc_arr(i,j,k).y=net.param.r_dc/2;
            else
                r_dc_arr(i,j,k).y=r_img(i,j,k);
            end
        end
    end
end

% Depth Resistors (z)
for i=1:Ni+1
    for j=1:Nj+1        
        for k=1:Nk
            if net.param.small_edges==1 && (k==1||k==Nk)
                r_dc_arr(i,j,k).z=net.param.r_dc/2;
            else
                r_dc_arr(i,j,k).z=r_img(i,j,k);
            end
        end
    end
end




