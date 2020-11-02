function plot_velec_fn(net,opt)

idx=[net.s.n1j];
idy=[net.s.n1i];
idz=[net.s.n1k];
idx2=[net.s.n2j];
idy2=[net.s.n2i];
idz2=[net.s.n2k];

if opt==1
    fontsz=13;
else
    fontsz=9;
end

s_arr=[net.s];
Ve_arr=abs([s_arr.Ve]);
Ve_max=max(Ve_arr,[],'all');

for i=1:length(net.s)
    idxm=idx(i);
    idzm=idz(i);
    if idxm==1
        idxm=idxm;
        %idym=idym-round(net.param.EIT_spacing/2);
    elseif idxm==net.Nj+1
        idxm=idx(i)-2.5;
        %idym=idym+round(net.param.EIT_spacing/2);
    elseif idzm==1
        idzm=idz(i)+0.5;
        idxm=idx(i)-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    elseif idzm==net.Nk+1
        idzm=idzm-0.5;
        %idxm=idxm+round(net.param.EIT_spacing/2);
    end
    if Ve_max<10
        text(idxm,idzm,strcat('Ve_{',num2str(i),'}=',num2str(abs(net.s(i).Ve)*1e3,'%1.0f')),'color','black','FontWeight','bold','FontSize',fontsz)
    else
        text(idxm,idzm,strcat('Ve_{',num2str(i),'}=',num2str(abs(net.s(i).Ve),'%1.1f')),'color','black','FontWeight','bold','FontSize',fontsz)
    end
end

