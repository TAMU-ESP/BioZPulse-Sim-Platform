function [net]=pulse_art_fn(net)


  
Nelec=net.param.Nelec;
EIT_spacing=net.param.EIT_spacing;
direction=net.param.direction;
i_dia_x=net.param.dia;
i_dia_y=net.param.dia;
s_dia_x=net.param.dia;
s_dia_y=net.param.dia;

Ni=20;  %Even
Nj=Nelec/4*EIT_spacing;
Nk=Nelec/4*EIT_spacing;
% Create Electrode Array
for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk+1
            node_arr(i,j,k)=xyz2node_fn(i,j,k,Ni,Nj,Nk);
        end
    end
end

S=size(r_var_freq_arr,2);
r_art_dc=4/7*r_dc;
ri_art_dc=4/7*ri_dc;
ci_art_dc=1/(4/7)*ci_dc;
r_var_freq_arr=[1];
r_var_time_arr=1./min(r_var_freq_arr,[],2)*2;


r_var_delay_arr=50e-3*ones(1,S);
r_var_n1j=[round((Nj)*(3/4))];    % R var start location vertical
r_var_n1i=[round((Ni)*(3/4))];    % R var start location horizontal
r_var_n1k_arr=[1]*ones(1,S);     % R var start location depth
r_var_max=0;
ri_var_max=ri_art_dc*1e-2*ones(1,S);
ci_var_max=ci_art_dc*1e-2*ones(1,S);
r_var_steps=3*ones(1,S);
r_var_v=0*ones(1,S);       % 1-> vertical, 0-> horizontal
r_var_dim=2*ones(1,S);     % Dimension (1,2,3)
r_var_num=Ni*ones(1,S);     % Number of resistors (1,..)
r_var_nan=1*ones(1,S);     % 0-> No NaN, 1-> NaN
r_var_seg=Ni*ones(1,S);

eit_elec_array.i_src_n1i_arr=i_src_n1i_arr;
eit_elec_array.i_src_n1j_arr=i_src_n1j_arr;
eit_elec_array.i_src_n1k_arr=i_src_n1k_arr;
eit_elec_array.i_src_n2i_arr=i_src_n2i_arr;
eit_elec_array.i_src_n2j_arr=i_src_n2j_arr;
eit_elec_array.i_src_n2k_arr=i_src_n2k_arr;

eit_elec_array.s_n1i_arr=s_n1i_arr;
eit_elec_array.s_n1j_arr=s_n1j_arr;
eit_elec_array.s_n1k_arr=s_n1k_arr;
eit_elec_array.s_n2i_arr=s_n2i_arr;
eit_elec_array.s_n2j_arr=s_n2j_arr;
eit_elec_array.s_n2k_arr=s_n2k_arr;
net.Ni=Ni;
net.Nj=Nj;
net.Nk=Nk;

net.eit_elec_array=eit_elec_array;

