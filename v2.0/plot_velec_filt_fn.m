function plot_velec_filt_fn(net)

idx=[net.s.n1j];
idy=[net.s.n1i];

for i=1:length(net.s)
    idxm=idx(i);
    idym=idy(i);
    if idxm==1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idym=idym-round(net.param.EIT_spacing/2);
    elseif idxm==net.Nj+1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idym=idym+round(net.param.EIT_spacing/2);
    elseif idym==1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    elseif idym==net.Ni+1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    end
    plot(idxm,idym,'ro','markersize',8,'linewidth',2);
end

idx=[net.s.n2j];
idy=[net.s.n2i];

for i=1:length(net.s)
    idxm=idx(i);
    idym=idy(i);
    if idxm==1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idym=idym-round(net.param.EIT_spacing/2);
    elseif idxm==net.Nj+1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idym=idym+round(net.param.EIT_spacing/2);
    elseif idym==1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    elseif idym==net.Ni+1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    end
    plot(idxm,idym,'ro','markersize',8,'linewidth',2);
end

idx=[net.curr.n1j];
idy=[net.curr.n1i];

for i=1:length(net.curr)
    idxm=idx(i);
    idym=idy(i);
    if idxm==1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idym=idym-round(net.param.EIT_spacing/2);
    elseif idxm==net.Nj+1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idym=idym+round(net.param.EIT_spacing/2);
    elseif idym==1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    elseif idym==net.Ni+1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    end
    plot(idxm,idym,'go','markersize',8,'linewidth',2);
end

idx=[net.curr.n2j];
idy=[net.curr.n2i];

for i=1:length(net.curr)
    idxm=idx(i);
    idym=idy(i);
    if idxm==1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idym=idym-round(net.param.EIT_spacing/2);
    elseif idxm==net.Nj+1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idym=idym+round(net.param.EIT_spacing/2);
    elseif idym==1
        idym=idym-0.5;
        idxm=idxm-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    elseif idym==net.Ni+1
        idxm=idxm-0.5;
        idym=idym-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    end
    plot(idxm,idym,'go','markersize',8,'linewidth',2);
end