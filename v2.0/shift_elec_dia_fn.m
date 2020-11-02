function [nim,njm]=shift_elec_dia_fn(ni,nj,Ni,Nj,dia_x,dia_y)

nim=ni;
njm=nj;
if ni==0
    njm=nj-round(dia_x/2);
elseif nj==Nj
    nim=ni-round(dia_y/2);
    njm=nj-dia_x;
elseif ni==Ni
    nim=ni-dia_y;
    njm=nj-round(dia_x/2);
elseif nj==0
    nim=ni-round(dia_y/2);
end