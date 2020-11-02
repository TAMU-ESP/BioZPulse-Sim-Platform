function [net]=eit_elec_array_fn(net,i_num_in)



Nelec=net.param.Nelec;
EIT_spacing=net.param.EIT_spacing;
direction=net.param.direction;
i_dia_x=net.param.dia;
i_dia_y=net.param.dia;
s_dia_x=net.param.dia;
s_dia_y=net.param.dia;

Nj=Nelec/4*EIT_spacing;
Nk=Nelec/4*EIT_spacing;
Ni=net.Ni;
% Create Electrode Array
for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk+1
            node_arr(i,j,k)=xyz2node_fn(i,j,k,Ni,Nj,Nk);
        end
    end
end
e=EIT_spacing;
if direction==0
    %i_src_n1_arr_all_ini=[squeeze(node_arr(round(Nk/2)+1,1,e/2+1:e:end-(e/2)))',node_arr(round(Nk/2)+1,e/2+1:e:end-(e/2),end),squeeze(node_arr(round(Nk/2)+1,end,end-(e/2):-e:e/2+1))',node_arr(round(Nk/2)+1,end-e/2:-e:e/2+1,1)]';
    i_src_n1_arr_all_ini=[node_arr(round(Ni/2)+1,e/2+1:e:end-(e/2),1),squeeze(node_arr(round(Ni/2)+1,end,e/2+1:e:end-(e/2)))',node_arr(round(Ni/2)+1,end-e/2:-e:e/2+1,end),squeeze(node_arr(round(Ni/2)+1,1,end-(e/2):-e:e/2+1))']';
    
    %i_src_n1_arr_all_ini2=fliplr(i_src_n1_arr_all_ini);
    %i_src_n1_arr_all=circshift(i_src_n1_arr_all_ini2,-1*round(length(i_src_n1_arr_all_ini2)/4));
    i_src_n1_arr_all=(i_src_n1_arr_all_ini);
    i_src_n2_arr_all=circshift(i_src_n1_arr_all,-1);
    i_src_n1_arr=i_src_n1_arr_all(1:end);
    i_src_n2_arr=i_src_n2_arr_all(1:end);
    if ~exist('i_num_in')
        i_num_in=length(i_src_n1_arr);
    end
    for t1=1:min(length(i_src_n1_arr),i_num_in)
        i_src_n1=i_src_n1_arr(t1);
        i_src_n2=i_src_n2_arr(t1);
        [i_src_n1i,i_src_n1j,i_src_n1k]=ind2sub(size(node_arr),find((node_arr)==i_src_n1));
        [i_src_n2i,i_src_n2j,i_src_n2k]=ind2sub(size(node_arr),find((node_arr)==i_src_n2));
        i_src_n1i=i_src_n1i;
        i_src_n2i=i_src_n2i;
        i_src_n1j=i_src_n1j;
        i_src_n2j=i_src_n2j;
        i_src_n1k=i_src_n1k;
        i_src_n2k=i_src_n2k;
        
        %s_all=[squeeze(node_arr(1,round(e/2)+1:e:end-round(e/2),:)),squeeze(permute(node_arr(round(e/2)+1:e:end-round(e/2),end,:),[2 1 3 4])),squeeze(node_arr(end,end-round(e/2):-e:round(e/2)+1,:)),squeeze(permute(node_arr(end-round(e/2):-e:round(e/2)+1,1,:),[2 1 3 4]))]-1;
        %         s_all_no_i=s_all(s_all~=i_src_n1);
        %         s_all_no_i=s_all_no_i(s_all_no_i~=i_src_n2);
        %         s_all2=circshift(s_all,-(t1+1));
        %         s_n1=s_all2(1:length(i_src_n1_arr)-3);
        %         s_n2=s_all2(2:length(i_src_n1_arr)-3+1) ;
        %         s_n1=circshift(s_n1,+(t1+1-2));
        %         s_n2=circshift(s_n2,+(t1+1-2));
        
        s_all=i_src_n1_arr_all;
        s_pair_valid=s_all;
        s_pair_valid((s_all==i_src_n1)|(s_all==i_src_n2))=NaN;
        s_pair_valid_arr(1,:)=s_pair_valid;
        s_pair_valid_arr(2,:)=circshift(s_pair_valid,-1);
        s_pair_valid_arr2=isnan(s_pair_valid_arr);
        s_pair_valid_arr3=s_pair_valid_arr2(1,:)|s_pair_valid_arr2(2,:);
        s_n1=s_pair_valid_arr(1,:);
        s_n1(s_pair_valid_arr3)=[];
        s_n2=s_pair_valid_arr(2,:);
        s_n2(s_pair_valid_arr3)=[];
        
        for s=1:length(s_n1)
            [s_n1i(s),s_n1j(s),s_n1k(s)]=ind2sub(size(node_arr),find((node_arr)==s_n1(s)));
            [s_n2i(s),s_n2j(s),s_n2k(s)]=ind2sub(size(node_arr),find((node_arr)==s_n2(s)));
        end
        s_n1i=s_n1i;
        s_n2i=s_n2i;
        s_n1j=s_n1j;
        s_n2j=s_n2j;
        s_n1k=s_n1k;
        s_n2k=s_n2k;
        
        for m=1:length(i_src_n1i)
            [i_src_n1im(m),i_src_n1jm(m)]=shift_elec_dia_fn(i_src_n1i(m),i_src_n1j(m),Ni,Nj,i_dia_x,i_dia_y);
            [i_src_n2im(m),i_src_n2jm(m)]=shift_elec_dia_fn(i_src_n2i(m),i_src_n2j(m),Ni,Nj,i_dia_x,i_dia_y);
        end
        for m=1:length(s_n1i)
            [s_n1im(m),s_n1jm(m)]=shift_elec_dia_fn(s_n1i(m),s_n1j(m),Ni,Nj,s_dia_x,s_dia_y);
            [s_n2im(m),s_n2jm(m)]=shift_elec_dia_fn(s_n2i(m),s_n2j(m),Ni,Nj,s_dia_x,s_dia_y);
        end
        
        i_src_n1i_arr(t1,:)=i_src_n1im;
        i_src_n1j_arr(t1,:)=i_src_n1jm;
        i_src_n1k_arr(t1,:)=i_src_n1k;
        i_src_n2i_arr(t1,:)=i_src_n2im;
        i_src_n2j_arr(t1,:)=i_src_n2jm;
        i_src_n2k_arr(t1,:)=i_src_n2k;
        i_src_n1_arr(t1,:)=i_src_n1;
        i_src_n2_arr(t1,:)=i_src_n2;
        
        s_n1i_arr(t1,:)=s_n1im;
        s_n1j_arr(t1,:)=s_n1jm;
        s_n1k_arr(t1,:)=s_n1k;
        s_n2i_arr(t1,:)=s_n2im;
        s_n2j_arr(t1,:)=s_n2jm;
        s_n2k_arr(t1,:)=s_n2k;
        s_n1_arr(t1,:)=s_n1;
        s_n2_arr(t1,:)=s_n2;
    end
else
    i_src_n1_arr_all=[node_arr(round(Ni/2)+1,e/2+1:e:end-(e/2),1),squeeze(node_arr(round(Ni/2)+1,end,e/2+1:e:end-(e/2)))',node_arr(round(Ni/2)+1,end-e/2:-e:e/2+1,end),squeeze(node_arr(round(Ni/2)+1,1,end-(e/2):-e:e/2+1))']';
    i_src_n2_arr_all=circshift(i_src_n1_arr_all,length(i_src_n1_arr_all)/2);
    i_src_n1_arr=i_src_n1_arr_all(1:end);
    i_src_n2_arr=i_src_n2_arr_all(1:end);
    if ~exist('i_num_in')
        i_num_in=length(i_src_n1_arr);
    end
    for t1=1:min(length(i_src_n1_arr),i_num_in)
        i_src_n1=i_src_n1_arr(t1);
        i_src_n2=i_src_n2_arr(t1);
        [i_src_n1i,i_src_n1j,i_src_n1k]=ind2sub(size(node_arr),find((node_arr)==i_src_n1));
        [i_src_n2i,i_src_n2j,i_src_n2k]=ind2sub(size(node_arr),find((node_arr)==i_src_n2));
        i_src_n1i=i_src_n1i;
        i_src_n2i=i_src_n2i;
        i_src_n1j=i_src_n1j;
        i_src_n2j=i_src_n2j;
        i_src_n1k=i_src_n1k;
        i_src_n2k=i_src_n2k;
        
        s_all=i_src_n1_arr_all;
        %s_all_no_i=s_all(s_all~=i_src_n1);
        %s_all_no_i=s_all_no_i(s_all_no_i~=i_src_n2);
        s_all2=circshift(s_all,-(t1));
        %s_n1=s_all2(1:length(i_src_n1_arr_all)-2);
        s_n2=repmat(s_all2(1),1,length(i_src_n1_arr_all)-3);
        s_n2=s_n2(s_n2~=i_src_n2);
        s_n1=s_all2(2:length(i_src_n1_arr_all)-2+1);
        s_n1=s_n1(s_n1~=i_src_n2);
        
        for s=1:length(s_n1)
            [s_n1i(s),s_n1j(s),s_n1k(s)]=ind2sub(size(node_arr),find((node_arr)==s_n1(s)));
            [s_n2i(s),s_n2j(s),s_n2k(s)]=ind2sub(size(node_arr),find((node_arr)==s_n2(s)));
        end
        s_n1i=s_n1i;
        s_n2i=s_n2i;
        s_n1j=s_n1j;
        s_n2j=s_n2j;
        s_n1k=s_n1k;
        s_n2k=s_n2k;
        
        for m=1:length(i_src_n1i)
            [i_src_n1im(m),i_src_n1jm(m)]=shift_elec_dia_fn(i_src_n1i(m),i_src_n1j(m),Ni,Nj,i_dia_x,i_dia_y);
            [i_src_n2im(m),i_src_n2jm(m)]=shift_elec_dia_fn(i_src_n2i(m),i_src_n2j(m),Ni,Nj,i_dia_x,i_dia_y);
        end
        for m=1:length(s_n1i)
            [s_n1im(m),s_n1jm(m)]=shift_elec_dia_fn(s_n1i(m),s_n1j(m),Ni,Nj,s_dia_x,s_dia_y);
            [s_n2im(m),s_n2jm(m)]=shift_elec_dia_fn(s_n2i(m),s_n2j(m),Ni,Nj,s_dia_x,s_dia_y);
        end
        
        i_src_n1i_arr(t1,:)=i_src_n1im;
        i_src_n1j_arr(t1,:)=i_src_n1jm;
        i_src_n1k_arr(t1,:)=i_src_n1k;
        i_src_n2i_arr(t1,:)=i_src_n2im;
        i_src_n2j_arr(t1,:)=i_src_n2jm;
        i_src_n2k_arr(t1,:)=i_src_n2k;
        i_src_n1_arr(t1,:)=i_src_n1;
        i_src_n2_arr(t1,:)=i_src_n2;
        
        s_n1i_arr(t1,:)=s_n1im;
        s_n1j_arr(t1,:)=s_n1jm;
        s_n1k_arr(t1,:)=s_n1k;
        s_n2i_arr(t1,:)=s_n2im;
        s_n2j_arr(t1,:)=s_n2jm;
        s_n2k_arr(t1,:)=s_n2k;
        s_n1_arr(t1,:)=s_n1;
        s_n2_arr(t1,:)=s_n2;
    end
end

elec_array.i_src_n1i_arr=i_src_n1i_arr;
elec_array.i_src_n1j_arr=i_src_n1j_arr;
elec_array.i_src_n1k_arr=i_src_n1k_arr;
elec_array.i_src_n2i_arr=i_src_n2i_arr;
elec_array.i_src_n2j_arr=i_src_n2j_arr;
elec_array.i_src_n2k_arr=i_src_n2k_arr;

elec_array.s_n1i_arr=s_n1i_arr;
elec_array.s_n1j_arr=s_n1j_arr;
elec_array.s_n1k_arr=s_n1k_arr;
elec_array.s_n2i_arr=s_n2i_arr;
elec_array.s_n2j_arr=s_n2j_arr;
elec_array.s_n2k_arr=s_n2k_arr;
elec_array.i_src_n1_arr=i_src_n1_arr;
elec_array.i_src_n2_arr=i_src_n2_arr;
elec_array.s_n1_arr=s_n1_arr;
elec_array.s_n2_arr=s_n2_arr;
net.Ni=Ni;
net.Nj=Nj;
net.Nk=Nk;

net.elec_array=elec_array;
net.param.i_dia_x=i_dia_x;
net.param.s_dia_x=s_dia_x;
net.param.i_dia_y=i_dia_y;
net.param.s_dia_y=s_dia_y;