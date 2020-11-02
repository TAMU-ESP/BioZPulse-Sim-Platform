function net_plot_fn(net_arr)

n=size(net_arr,1);
m=size(net_arr,2);

if (n==1 && m==1) || (net_arr{1,1}(1).param.plot_eit==2 && net_arr{1,1}(1).param.debug==0)
    for i=1:n
        for j=1:m
            V_node_all=[net_arr{i,j}.V_node];
            V_node_max=max(max(V_node_all));
            V_node_min=min(min(V_node_all));
            for k=1:length(net_arr{i})
                net_arr{i,j}(k).param.plot_path=strcat(net_arr{i,j}(k).param.plot_path,'/',num2str(i),',',num2str(i),'_');
                spice_netlist_eit_plot_fn(net_arr{i,j}(k),V_node_min,V_node_max);
            end
        end
    end
end