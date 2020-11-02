function [dVe_arr_all,dVenorm_arr_all]=sens_calc_fn(net_arr,net_ref)

curr_idx_arr=1:length(net_arr{1,1});
s_idx_arr=1:length(net_arr{1,1}(1).s);
for k=1:length(curr_idx_arr)
    for s=1:length(s_idx_arr)
        curr_idx=curr_idx_arr(k);
        s_idx=s_idx_arr(s);
        n=size(net_arr,1);
        m=size(net_arr,2);
        for i=1:n
            for j=1:m
                s_arr=cat(1,net_arr{i,j}(curr_idx).s);
                dVe_arr_all(k,s,j,i)=s_arr(s_idx).Vs-net_ref{1}(curr_idx).s(s_idx).Vs;
            end
        end
    end
end

% Combined senestivity map
dVe_arr_all_re=reshape(dVe_arr_all,[],n,m);
dVenorm_arr_all=squeeze(vecnorm(dVe_arr_all_re,2,1));