function obj_fn = myMPC_ObjFn(Uf_vec)

global N_pred N_con We Wdu n_st n_op n_ip
global xhat_k ek_f phy gama Lp_inf C_mat rk_pred uk_minus_1 ek_filt

zk_pred = zeros(n_st, N_pred);

zk_pred(:,1) = xhat_k;
yk_pred = zeros(n_op, N_pred);
yk_pred(:,1) = C_mat*zk_pred(:,1) + ek_filt;

ukf = zeros(n_ip, N_pred);
n2 = 0;
for k=1:N_pred
    if (k <= N_con)
        n1 = n2 + 1; n2 = k*n_ip; % n1 increase by 1 and n2 increases by multiple of 2 
        ukf(:,k) = Uf_vec(n1:n2,1)';
    else
        ukf(:,k) = Uf_vec(n1:n2,1)';
    end
end

del_uk = ukf(:,1) - uk_minus_1;
err_pred = rk_pred - yk_pred(:,1);

obj_fn = err_pred'*We*err_pred + del_uk'*Wdu*del_uk;
for k =2:N_pred
    zk_pred(:,k) = phy*zk_pred(:,k-1)+gama*ukf(:,k-1)+Lp_inf*ek_filt;
    yk_pred(:,k) = C_mat*zk_pred(:,k) + ek_filt;
    err_pred = rk_pred - yk_pred(:,k);
    del_uk = ukf(:,k) - ukf(:,k-1);
    obj_fn = obj_fn + err_pred'*We*err_pred+del_uk'*Wdu*del_uk;
end