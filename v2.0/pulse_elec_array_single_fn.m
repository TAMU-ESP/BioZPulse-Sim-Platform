function [net]=pulse_elec_array_single_fn(net)

Ni=(net.Ni);Nj=(net.Nj);Nk=(net.Nk);
dia_x=round(net.param.Ex/net.param.L);
dia_y=round(net.param.Ey/net.param.L);
pi_x=round(net.param.PIx/net.param.L);
pi_y=round(net.param.PIy/net.param.L);
si_y=round(net.param.SIy/net.param.L);
pv_x=round(net.param.PVx/net.param.L);
pv_y=round(net.param.PVy/net.param.L);
sv_y=round(net.param.SVy/net.param.L);
s_num=length(sv_y);


i_src_n1i=pi_y+1-round(0.5*dia_y);
i_src_n1j=pi_x+1-round(dia_x/2);
i_src_n2i=i_src_n1i+si_y;
i_src_n2j=i_src_n1j;
i_src_n1k=1;
i_src_n2k=1;

s_n1i=[];s_n1j=[];s_n2i=[];s_n2j=[];
if s_num>0
    for i=1:s_num        
        s_n1i(i)=pv_y(i)+1-round(0.5*dia_y);
        s_n2i(i)=s_n1i(i)+sv_y(i);
        s_n1j(i)=pv_x(i)+1-round(dia_x/2);
        s_n2j(i)=s_n1j(i);
        s_n1k(i)=1;
        s_n2k(i)=1;
    end
end


elec_array.i_src_n1i_arr=i_src_n1i;
elec_array.i_src_n1j_arr=i_src_n1j;
elec_array.i_src_n1k_arr=i_src_n1k;
elec_array.i_src_n2i_arr=i_src_n2i;
elec_array.i_src_n2j_arr=i_src_n2j;
elec_array.i_src_n2k_arr=i_src_n2k;

elec_array.s_n1i_arr=s_n1i;
elec_array.s_n1j_arr=s_n1j;
elec_array.s_n1k_arr=s_n1k;
elec_array.s_n2i_arr=s_n2i;
elec_array.s_n2j_arr=s_n2j;
elec_array.s_n2k_arr=s_n2k;

net.elec_array=elec_array;
net.param.i_dia_x=dia_x;
net.param.s_dia_x=dia_x;
net.param.i_dia_y=dia_y;
net.param.s_dia_y=dia_y;

