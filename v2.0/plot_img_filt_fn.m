function plot_img_filt_fn(eit,elec_arr,opt)
% opt=0 -> No electrode, 1-> with electrodes , 2-> with electrodes & arrows

n=size(eit,1);
m=size(eit,2);

%eit{1,1}(1).param.plot_eit==1 &&
if (n*m<16) || (eit{1,1}.net(1).param.debug==1)
    for i=1:n
        for j=1:m
            r_img_filt_all(i,j,:,:)=[eit{i,j}.net(1).r_img_filt];
        end
    end
    clims=[min(min(min(min(r_img_filt_all)))) max(max(max(max(r_img_filt_all))))];
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
                    mesh_plot_fn(squeeze(eit{i,j}.net(p).r_img_filt),1,clims);
                    axis equal;
                    xlim([0.5 eit{i,j}.net(1).Nj+1.5-1])
                    ylim([0.5 eit{i,j}.net(1).Ni+1.5-1])
                elseif opt==1
                    mesh_plot_fn(squeeze(eit{i,j}.net(p).r_img_filt(:,:,1)),1,clims);
                    axis equal;
                    plot_velec_filt_fn(eit{i,j}.net(p))
                    xlim([0.5 eit{i,j}.net(1).Nj+1.5-1])
                    ylim([0.5 eit{i,j}.net(1).Ni+1.5-1])
                else
                    mesh_plot_fn(squeeze(eit{i,j}.net(p).img_filt_arr(1:end-1,1:end-1,1)),1,clims);
                    axis equal;
                    Vs_arrow_plot_fn(eit{i,j}.net(p));
                    xlim([0.5 eit{i,j}.net(1).Nj+1.5-1])
                    ylim([0.5 eit{i,j}.net(1).Ni+1.5-1])
                end
                title(strcat('Orig. Conductivity Image ',num2str(i),',',num2str(j),' (S)'));
                xlabel('X');ylabel('Y');              
                ax = gca; ax.FontSize=11;
            end
        end
    end
    
    if eit{1,1}.net(1).param.write_figures==1 && (eit{1,1}.net(1).param.debug==0||eit{1,1}.net(1).param.eidors_calc_en==1)
        if ~exist(eit{1,1}.net(1).param.path,'dir')
            mkdir(eit{1,1}.net(1).param.path);
        end
        file_str=strcat(eit{1,1}.net(1).param.path);
        print(strcat(file_str,'img_eidors_input',num2str(n),'x',num2str(m),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'img_eidors_input',num2str(n),'x',num2str(m),'.fig')); %<-Save as PNG with 300 DP
    end
end


if (n*m<16) || (eit{1,1}.net(1).param.eidors_calc_en==1)
    for i=1:n
        for j=1:m
            r_img_reconst_all(i,j,:,:)=[eit{i,j}.inv.img_reconst];
        end
    end
    clims=[min(min(min(min(r_img_reconst_all)))) max(max(max(max(r_img_reconst_all))))];
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
            [c d]=size(eit{i,j}.inv.img_reconst);
            %for k=1:length(elec_arr)
                p=elec_arr(k);
                ind=ind+1;
                subplot(a,b,ind)
                if opt==0
                    mesh_plot_fn(squeeze(eit{i,j}.inv.img_reconst),1,clims);
                    axis equal;
                    xlim([0.5 d+1.5-1])
                    ylim([0.5 c+1.5-1])
                elseif opt==1
                    mesh_plot_fn(squeeze(eit{i,j}.inv.img_reconst),1,clims);
                    axis equal;
                    plot_velec_filt_fn(eit{i,j}.net(p))
                    xlim([0.5 d+1.5-1])
                    ylim([0.5 c+1.5-1])
                else
                    mesh_plot_fn(squeeze(eit{i,j}.inv.img_reconst(1:end-1,1:end-1,1)),1,clims);
                    axis equal;
                    Vs_arrow_plot_fn(eit{i,j}.net(p));
                    xlim([0.5 d+1.5-1])
                    ylim([0.5 c+1.5-1])
                end
                title(strcat('Reconst. Conductivity Image ',num2str(i),',',num2str(j),' (S)'));
                xlabel('X');ylabel('Y');              
                ax = gca; ax.FontSize=11;
            %end
        end
    end
    
    if eit{1,1}.net(1).param.write_figures==1 && (eit{1,1}.net(1).param.debug==0||eit{1,1}.net(1).param.eidors_calc_en==1)
        if ~exist(eit{1,1}.net(1).param.path,'dir')
            mkdir(eit{1,1}.net(1).param.path);
        end
        file_str=strcat(eit{1,1}.net(1).param.path);
        print(strcat(file_str,'img_eidors_reconst',num2str(n),'x',num2str(m),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'img_eidors_reconst',num2str(n),'x',num2str(m),'.fig')); %<-Save as PNG with 300 DP
    end
end