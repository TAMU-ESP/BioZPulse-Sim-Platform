function plot_img_arr_fn(net_arr,opt,clims_in,colormap_in)
% opt=1 => V,I   opt=2 => V,I,Z
n=size(net_arr,1);
m=size(net_arr,2);

for i=1:n
    for j=1:m
        %for k=1:length(elec_arr)
        %for y=1:length(freq_arr)
        %p=elec_arr(k);
        %f=freq_arr(y);
        %net_arr{i,j}(p,f).param.plot_path=strcat(net_arr{i,j}(p,f).param.path);
        %net_arr{i,j}(p,f).plot.freq_arr=freq_val_arr(freq_arr);
        %net_arr{i,j}(p,f).plot.freq_ind=y;
        %for t=y_arr
        plot_single_img_fn(net_arr{i,j}(1),opt,clims_in,colormap_in)
        if net_arr{i,j}(1).param.write_figures==1 && net_arr{i,j}(1).param.debug==0            
            file_str=strcat(net_arr{i,j}(1).param.path);
            pause(1);
            print(strcat(file_str,'spice_img_',num2str(i),'_',num2str(j),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
            %savefig(strcat(file_str,'spice_2D_',num2str(i),'_',num2str(j),'_',num2str(p),'_',num2str(f),'_',num2str(y),'.fig')); %<-Save as PNG with 300 DP
        end
        %end
        %end
        % end
    end
end


