function net=img_update_fn(net)

net.ri_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)=net.ri_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)*(1+net.param.update.delta_r_dc);  % percentage
net.re_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)=net.re_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)*(1+net.param.update.delta_r_dc);  % percentage
net.cm_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)=net.cm_img(net.param.update.center_y,net.param.update.center_x,net.param.update.center_z)*(1-net.param.update.delta_r_dc);  % percentage
