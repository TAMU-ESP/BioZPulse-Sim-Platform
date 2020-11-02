function [net]=pulse_elec_array_fn(net)

Nk=10;
dia_arr=[5];
dia_ar=1;
s_num=1;
margin_x=15;
margin_y=6;
spacing_x=max(dia_arr)*dia_ar+8;
spacing_y=max(dia_arr)+[3];
%spacing_y_sns_arr=max(dia_arr)+[8 16 24];
%spacing_y_sns_arr=max(dia_arr)+[10 20 30];
spacing_y_sns_arr=max(dia_arr)+[3];
dia_x=dia_arr*dia_ar;
dia_y=dia_arr;

Nj=(max(dia_arr)*dia_ar+2*margin_x);
i_src_n1j_arr=round(Nj*0.25-dia_x/2);

t1=1;
spacing_y_sns=spacing_y_sns_arr(t1);
Ni=(max(dia_arr)+(s_num*2-1)*spacing_y_sns+2*spacing_y+2*margin_y);
spacing_arr_in=spacing_y;
i_src_n1j=i_src_n1j_arr;
% Sensing

i_src_n1i=round((Ni-(dia_y+(s_num*2-1)*spacing_y_sns+2*spacing_y))/2)+1;
s_n1i=[];s_n1j=[];s_n2i=[];s_n2j=[];
if s_num>0
    for i=1:s_num
        if s_num==1
            s_n1i(i)=i_src_n1i+spacing_y*(2*i-1);
            s_n2i(i)=s_n1i(i)+spacing_y_sns;
            s_n1j(i)=i_src_n1j;
            s_n2j(i)=s_n1j(i);
            s_n1k(i)=1;
            s_n2k(i)=1;
        else
            s_n1i(i)=s_n2i(i-1)+spacing_y*(2*i-1);
            s_n2i(i)=s_n1i(i)+spacing_y_sns;
            s_n1j(i)=i_src_n1j;
            s_n2j(i)=s_n1j(i);
            s_n1k(i)=1;
            s_n2k(i)=1;
        end
    end
end
i_src_n2i=s_n2i(end)+spacing_y;
i_src_n2j=i_src_n1j;
i_src_n1k=1;
i_src_n2k=1;

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

net.Ni=Ni;
net.Nj=Nj;
net.Nk=Nk;
net.elec_array=elec_array;
net.param.i_dia_x=dia_x;
net.param.s_dia_x=dia_x;
net.param.i_dia_y=dia_y;
net.param.s_dia_y=dia_y;

