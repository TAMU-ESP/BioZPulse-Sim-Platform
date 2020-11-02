function plot_spice_arr_fn(net_arr,elec_arr,freq_arr,y_arr,opt)
% opt=1 => V,I   opt=2 => V,I,Z
n=size(net_arr,1);
m=size(net_arr,2);

factor=0.75;
V_node_all=[];
for i=1:n
    for j=1:m
        V_node_all=[V_node_all;abs([net_arr{i,j}(elec_arr,freq_arr).V_node])];
    end
end
V_node_max=max(V_node_all,[],'all')*factor;
%V_node_min=min(V_node_all,[],'all');
V_node_min=V_node_max*(1-factor);

Z_img_all=[];
for i=1:n
    for j=1:m
        Z_img_all=[Z_img_all;abs([net_arr{i,j}(elec_arr,freq_arr).Y_img])];
    end
end
%Z_img_all=abs([net_arr{:}(elec_arr,freq_arr).Y_img]);
Z_img_max=(max(Z_img_all,[],'all'));
Z_img_min=(min(Z_img_all,[],'all'));

ac_freq=net_arr{1,1}(1,1).param.ac_freq;
freq_val_arr=linspace(ac_freq(2),ac_freq(3),ac_freq(1));
for i=1:n
    for j=1:m
        for k=1:length(elec_arr)
            for y=1:length(freq_arr)
                p=elec_arr(k);
                f=freq_arr(y);
                %net_arr{i,j}(p,f).param.plot_path=strcat(net_arr{i,j}(p,f).param.path,num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_');
                net_arr{i,j}(p,f).param.plot_path=strcat(net_arr{i,j}(p,f).param.path);
                net_arr{i,j}(p,f).plot.freq_arr=freq_val_arr(freq_arr);
                net_arr{i,j}(p,f).plot.freq_ind=y;
                %plot_spice_fn(net_arr{i,j}(p,f),opt,V_node_min,V_node_max,Z_img_max,Z_img_min);
                for y=y_arr
                    plot_spice_v_i_2d_fn(net_arr{i,j}(p,f),opt,y,V_node_min,V_node_max);
                    if net_arr{i,j}(p,f).param.write_figures==1 && net_arr{i,j}(p,f).param.debug==0
                        file_str=strcat(net_arr{i,j}(p,f).param.plot_path);
                        print(strcat(file_str,'spice_2D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
                        savefig(strcat(file_str,'spice_2D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.fig')); %<-Save as PNG with 300 DP
                    end
                    plot_spice_v_i_3d_fn(net_arr{i,j}(p,f),opt,V_node_min,V_node_max);                         
                    if net_arr{i,j}(p,f).param.write_figures==1 && net_arr{i,j}(p,f).param.debug==0
                        file_str=strcat(net_arr{i,j}(p,f).param.plot_path);
                        print(strcat(file_str,'spice_3D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
                        savefig(strcat(file_str,'spice_3D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.fig')); %<-Save as PNG with 300 DP
                    end
                end
            end
        end
    end
end


