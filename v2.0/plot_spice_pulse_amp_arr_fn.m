function plot_spice_pulse_amp_arr_fn(net_arr,elec_arr,freq_arr,opt,factor)
% opt=1 => V,I   opt=2 => V,I,Z

n=size(net_arr,1);
m=size(net_arr,2);

for i=1:n
    for j=1:m
        K=size(net_arr{i,j}(1,1).V_node,1);
        V=size(net_arr{i,j}(1,1).V_node,3);
        for k=1:K
            for v=1:V
                [PAT(k,v),dVpp(k,v),Vdc(k,v)]=find_delay_fn(net_arr{i,j}(1,1).t_var,squeeze(abs(net_arr{i,j}(1,1).V_node(k,round(mean(net_arr{i,j}(1,1).curr.n1j)),v,:))),net_arr{i,j}(1,1).param.r_var_freq(1),0);
            end
        end
    end
end

V_node_all=[];
for i=1:n
    for j=1:m
        
        net_arr{i,j}(elec_arr,freq_arr).V_node=net_arr{i,j}(elec_arr,freq_arr).V_node(:,:,2:end,:);     
        net_arr{i,j}(elec_arr,freq_arr).Nk=net_arr{i,j}(elec_arr,freq_arr).Nk-1;
        net_arr{i,j}(elec_arr,freq_arr).img_arr=net_arr{i,j}(elec_arr,freq_arr).img_arr(:,:,2:end,:);
        if net_arr{i,j}(elec_arr,freq_arr).param.icalc_en==1
            net_arr{i,j}(elec_arr,freq_arr).i_vector=net_arr{i,j}(elec_arr,freq_arr).i_vector(:,:,2:end,:,:,:);
        end
        V_node_all=[V_node_all;abs([net_arr{i,j}(elec_arr,freq_arr).V_node(:,:,:,1)])];
    end
end
V_node_max1=max(dVpp,[],'all');
V_node_min1=min(dVpp,[],'all');
V_node_mean=mean([V_node_max1,V_node_min1]);

V_node_max=V_node_mean+(V_node_max1-V_node_mean)*factor;
V_node_min=V_node_mean-(V_node_mean-V_node_min1)*factor;

ac_freq=net_arr{1,1}(1,1).param.ac_freq;
freq_val_arr=linspace(ac_freq(2),ac_freq(3),ac_freq(1));
for i=1:n
    for j=1:m
        for k=1:length(elec_arr)
            for y=1:length(freq_arr)
                p=elec_arr(k);
                f=freq_arr(y);
                net_arr{i,j}(p,f).param.plot_path=strcat(net_arr{i,j}(p,f).param.path);
                net_arr{i,j}(p,f).plot.freq_arr=freq_val_arr(freq_arr);
                net_arr{i,j}(p,f).plot.freq_ind=y;
                %for t=y_arr
                plot_spice_pulse_amp_2d_fn(net_arr{i,j}(p,f),opt,round(mean(net_arr{i,j}(p,f).curr.n1j)));
                if net_arr{i,j}(p,f).param.write_figures==1 && net_arr{i,j}(p,f).param.debug==0
                    file_str=strcat(net_arr{i,j}(p,f).param.path);
                    print(strcat(file_str,'spice_2D_no_skin_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
                    %savefig(strcat(file_str,'spice_2D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.fig')); %<-Save as PNG with 300 DP
                end
                
                %end
            end
        end
    end
end


