function net=img_gen_fn(net)

%   Bassem Ibrahim
%   2020.11.02
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari
%
%   function r_dc_arr=imp_image_gen_fn(r_in,net)
%   Adjust bio-impoedance values to MATLAB

net.param.ltspice_path='C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
net.param.skin_elec_data_path='./data_out.mat';
net.param.cole_data_path='./cole_v4.mat';

if (net.param.debug==0)&&(net.param.write_figures==1 || isempty(strmatch(net.param.write_netlist,'0')))
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
else
    net.param.path='./';
end

load(net.param.cole_data_path);
cole_name={cole.name};
small_edges=0;

Ni=round(net.param.Ly/net.param.L);Nj=round(net.param.Lx/net.param.L);Nk=round(net.param.Lz/net.param.L);
net.param.small_edges=small_edges;

% Create R Image
ri_img=ones(Ni+1,Nj+1,Nk+1);
re_img=ones(Ni+1,Nj+1,Nk+1);
cm_img=ones(Ni+1,Nj+1,Nk+1);
a=0;
net.param.re_var_max=0;
net.param.r_var_freq=0;
if ~isempty(net.param.obj)
    for b=1:length(net.param.obj)
        if net.param.obj(b).art_en==1
            a=a+1;
            str_find=contains(cole_name,net.param.obj(b).tissue);
            net.param.ri_art_dc=cole(str_find).RI*net.param.obj(b).r_ratio;
            net.param.re_art_dc=cole(str_find).RE*net.param.obj(b).r_ratio;
            net.param.cm_art_dc=cole(str_find).CM/net.param.obj(b).r_ratio;
            net.param.r_var_freq(a)=[1];
            net.param.r_var_time(a)=1./min(net.param.r_var_freq(a),[],2)*2;
            net.param.r_var_delay(a)=14.7e-3;
            net.param.r_var_n1j(a)=floor((Nj)*net.param.obj(b).center.x-(Nj*net.param.obj(b).len.x/2));    % R var start location vertical
            net.param.r_var_n1i(a)=0;    % R var start location horizontal
            net.param.r_var_n1k(a)=floor((Nk)*net.param.obj(b).center.z-(Nk*net.param.obj(b).len.z/2));     % R var start location depth
            net.param.ri_var_max(a)=net.param.ri_art_dc*net.param.obj(b).r_art_ratio;
            net.param.re_var_max(a)=net.param.re_art_dc*net.param.obj(b).r_art_ratio;
            net.param.cm_var_max(a)=net.param.cm_art_dc*net.param.obj(b).r_art_ratio;
            net.param.r_var_steps(a)=3;
            net.param.r_var_v(a)=0;       % 1-> vertical, 0-> horizontal
            net.param.r_var_dim(a)=3;     % Dimension (1,2,3)
            net.param.r_var_num(a)=Ni;     % Number of resistors (1,..)
            net.param.r_var_nan(a)=1;     % 0-> No NaN, 1-> NaN
            net.param.r_var_seg(a)=Ni;
            net.param.r_var_dia(a)=max([floor(Nj*net.param.obj(b).len.x) 0]');
            net.param.r_var_tilt(a)=net.param.obj(b).tilt_en;
        else
            ix=floor(net.param.obj(b).center.x*(Nj)-net.param.obj(b).len.x/2*(Nj))+1:floor(net.param.obj(b).center.x*(Nj)+net.param.obj(b).len.x/2*(Nj))+1;
            iy=floor(net.param.obj(b).center.y*(Ni)-net.param.obj(b).len.y/2*(Ni))+1:floor(net.param.obj(b).center.y*(Ni)+net.param.obj(b).len.y/2*(Ni))+1;
            iz=floor(net.param.obj(b).center.z*(Nk)-net.param.obj(b).len.z/2*(Nk))+1:floor(net.param.obj(b).center.z*(Nk)+net.param.obj(b).len.z/2*(Nk))+1;
            %iz=round(net.param.obj(b).start.z*(Nk+1)):round((net.param.obj(b).start.z+net.param.obj(b).len.z)*(Nk+1));
            o(b).ix=ix;
            o(b).iy=iy;
            o(b).iz=iz;
            str_find=contains(cole_name,net.param.obj(b).tissue);
            if net.param.obj(b).tilt_en==0
                ri_img(iy(iy<(Ni+2)),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).RI*net.param.obj(b).r_ratio;
                re_img(iy(iy<(Ni+2)),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).RE*net.param.obj(b).r_ratio;
                cm_img(iy(iy<(Ni+2)),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).CM/net.param.obj(b).r_ratio;
            else
                ri_img(iy(iy<(Ni+2)&(iy>round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2))-net.param.obj(b).tilt_en)=cole(str_find).RI*net.param.obj(b).r_ratio;
                re_img(iy(iy<(Ni+2)&(iy>round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2))-net.param.obj(b).tilt_en)=cole(str_find).RE*net.param.obj(b).r_ratio;
                cm_img(iy(iy<(Ni+2)&(iy>round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2))-net.param.obj(b).tilt_en)=cole(str_find).CM/net.param.obj(b).r_ratio;
                
                ri_img(iy(iy<(Ni+2)&(iy<=round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).RI*net.param.obj(b).r_ratio;
                re_img(iy(iy<(Ni+2)&(iy<=round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).RE*net.param.obj(b).r_ratio;
                cm_img(iy(iy<(Ni+2)&(iy<=round((Ni+1)/2))),ix(0<ix&ix<(Nj+2)),iz(iz<(Nk+2)))=cole(str_find).CM/net.param.obj(b).r_ratio;
            end
        end
    end
    if ~isempty(net.param.obj)
        net.o=0;
    end
else
    net.o=0;
    net.param.re_var_max=0;
    net.param.r_var_freq=0;
end

net.ri_img=ri_img;
net.re_img=re_img;
net.cm_img=cm_img;

net.Ni=Ni;
net.Nj=Nj;
net.Nk=Nk;

