function plot_img_fn(net_arr,elec_arr,opt)
% opt=0 -> No electrode, 1-> with electrodes , 2-> with electrodes & arrows

n=size(net_arr,1);
m=size(net_arr,2);

%net_arr{1,1}(1).param.plot_eit==1 &&
if (n*m<16) || (net_arr{1,1}(1).param.debug==1)
    for i=1:n
        for j=1:m
            r_img_all(i,j,:,:,:)=[net_arr{i,j}(1).re_img];
        end
    end
    clims=[min(min(min(min(r_img_all)))) max(max(max(max(r_img_all))))];
    if clims(1)==clims(2)
        clims(1)=0;
        clims(2)=clims(2)*2;
    end
    figure
    a=ceil(sqrt(n*m*length(elec_arr)));
    b=ceil(n*m*length(elec_arr)/a);
    ind=0;
    for i=1:n
        for j=1:m
            for k=1:length(elec_arr)
                p=elec_arr(k);
                ind=ind+1;
                subplot(a,b,ind)
                if opt==0
                    mesh_plot_fn(squeeze(net_arr{i,j}(p).re_img),1,clims);
                elseif opt==1
                    mesh_plot_fn(squeeze(net_arr{i,j}(p).img_arr(round(net_arr{i,j}(p).Ni/2)+1,:,:)),1,clims);
                else
                    mesh_plot_fn(squeeze(net_arr{i,j}(p).img_arr(round(net_arr{i,j}(p).Ni/2)+1,:,:)),1,clims);
                    Vs_arrow_plot_fn(net_arr{i,j}(p));
                end
                title(strcat('Impedance Image ',num2str(i),',',num2str(j),',',num2str(p),' (\Omega)'));
                xlabel('X');ylabel('Z');
                axis equal;
                xlim([0.5 net_arr{i,j}(1).Nj+1.5])
                ylim([0.5 net_arr{i,j}(1).Nk+1.5])
                ax = gca; ax.FontSize=13;
            end
        end
    end
    
    if net_arr{1,1}(1).param.write_figures==1 && (net_arr{1,1}(1).param.debug==0)
        if ~exist(net_arr{1,1}(1).param.path,'dir')
            mkdir(net_arr{1,1}(1).param.path);
        end
        file_str=strcat(net_arr{1,1}(1).param.path);
        print(strcat(file_str,'img_all',num2str(n),'x',num2str(m),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'img_all',num2str(n),'x',num2str(m),'.fig')); %<-Save as PNG with 300 DP
    end
end


