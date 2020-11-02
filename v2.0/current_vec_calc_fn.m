function net=current_vec_calc_fn(net,opt)

% opt=0 => old, opt=1 => new

if net.param.icalc_en==1
    if opt==0
        
        V_node=net.V_node;
        % Current Distribution
        for i=1:net.Ni+1
            for j=1:net.Nj+1
                if i==1 && j==1
                    i1=0;
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=0;
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                elseif i==1 && j==net.Nj+1
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=0;
                    i3=0;
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                elseif j==1 && i==net.Ni+1
                    i1=0;
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=0;
                elseif j==net.Nj+1 && i==net.Ni+1
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=0;
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=0;
                elseif i==1
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=0;
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                elseif i==net.Ni+1
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=0;
                elseif j==1
                    i1=0;
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                elseif j==net.Nj+1
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=0;
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                else
                    i1=(V_node(i,j-1,1)-V_node(i,j,1))/net.grid.x.r_dc_arr(i,j-1,1);
                    i2=(V_node(i,j,1)-V_node(i,j+1,1))/net.grid.x.r_dc_arr(i,j,1);
                    i3=(V_node(i-1,j,1)-V_node(i,j,1))/net.grid.y.r_dc_arr(i-1,j,1);
                    i4=(V_node(i,j,1)-V_node(i+1,j,1))/net.grid.y.r_dc_arr(i,j,1);
                end
                i_vector(i,j,1,:)=[(i2+i1),(i3+i4)];
            end
        end
        
    else
        
        % Current Distribution
        for i=1:net.Ni+1
            for j=1:net.Nj+1
                for k=1:net.Nk+1
                    if i==1 && j==1 && k==1   %000
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==1 && j==net.Nj+1 && k==1 %010
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==1 && i==net.Ni+1 && k==1 %100
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==net.Nj+1 && i==net.Ni+1 && k==1 %110
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==1 && j==1 && k==net.Nk+1  %001
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif i==1 && j==net.Nj+1 && k==net.Nk+1  %011
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif j==1 && i==net.Ni+1 && k==net.Nk+1  %101
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif j==net.Nj+1 && i==net.Ni+1 && k==net.Nk+1  %111
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif i==1 && k==1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==net.Ni+1 && k==1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==1 && k==1
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==net.Nj+1 && k==1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==1 && k==net.Nk+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif i==net.Ni+1 && k==net.Nk+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif j==1 && k==net.Nk+1
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif j==net.Nj+1 && k==net.Nk+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    elseif i==1 && j==1
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==1 && j==net.Nj+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==net.Ni+1 && j==1
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==net.Ni+1 && j==net.Nj+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=0;
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==1
                        i1=0;
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif k==1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=0;
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif i==net.Ni+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=0;
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif j==net.Nj+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=0;
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    elseif k==net.Nk+1
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=0;
                    else
                        i1=net.I_node.x.Ri(i,j-1,k)+net.I_node.x.Re(i,j-1,k);
                        i2=net.I_node.x.Ri(i,j,k)+net.I_node.x.Re(i,j,k);
                        i3=net.I_node.y.Ri(i-1,j,k)+net.I_node.y.Re(i-1,j,k);
                        i4=net.I_node.y.Ri(i,j,k)+net.I_node.y.Re(i,j,k);
                        i5=net.I_node.z.Ri(i,j,k-1)+net.I_node.z.Re(i,j,k-1);
                        i6=net.I_node.z.Ri(i,j,k)+net.I_node.z.Re(i,j,k);
                    end
                    i_vector(i,j,k,:)=([(i2+i1),(i3+i4),(i5+i6)]);
                end
            end
        end
    end
    
    net.Z_node.x=(net.V_node(:,1:end-1,:,:,:)-net.V_node(:,2:end,:,:,:))./(net.I_node.x.Ri(:,1:end-1,:,:,:)+net.I_node.x.Re(:,1:end-1,:,:,:));
    net.Z_node.y=(net.V_node(1:end-1,:,:,:,:)-net.V_node(2:end,:,:,:,:))./(net.I_node.y.Ri(1:end-1,:,:,:,:)+net.I_node.y.Re(1:end-1,:,:,:,:));
    net.Z_node.z=(net.V_node(:,:,1:end-1,:,:)-net.V_node(:,:,2:end,:,:))./(net.I_node.z.Ri(:,:,1:end-1,:,:)+net.I_node.z.Re(:,:,1:end-1,:,:));
    if net.Ni>1
        net.Z_img=(net.Z_node.x(1:end-1,:,1:end-1,:,:)+net.Z_node.y(:,1:end-1,1:end-1,:,:)+net.Z_node.x(2:end,:,2:end,:,:)+net.Z_node.y(:,2:end,2:end,:,:)+net.Z_node.z(1:end-1,1:end-1,:,:,:)+net.Z_node.z(2:end,2:end,:,:,:))/6;
    else
        net.Z_img=(net.Z_node.x(1,:,1:end-1,:,:)+net.Z_node.x(1,:,2:end,:,:)+net.Z_node.z(1,1:end-1,:,:,:)+net.Z_node.z(1,2:end,:,:,:))/4;
    end
    net.Y_img=1./net.Z_img;
    
    net.i_vector=i_vector;
else
    net.i_vector=0;
end
