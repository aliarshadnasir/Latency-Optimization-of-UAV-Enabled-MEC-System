clc;clear;

dispi = 0;

op = 1000;

if (op == 100)
    Xmin = 0; Xmax = 100;
    Ymin = 0; Ymax = 100;
else
    Xmin = 0; Xmax = 1000;
    Ymin = 0; Ymax = 1000;
    csX = 1000; csY = 1000;
end
N_vec = [15:3:30];
result_vec = zeros(1,length(N_vec));
for N = N_vec
    
    M = 2;
    H = 150*ones(M,1);
    kap = 1e-27;
    g_tol = 1e-4;
    % w = zeros(N,2); % user position
    % w_30 = 100*rand(30,2)'
    % w_30 = [64.4318   37.8609   81.1580   53.2826   35.0727   93.9002   87.5943   55.0156   62.2475   58.7045   20.7742   30.1246   47.0923   23.0488   84.4309  19.4764   22.5922   17.0708   22.7664   43.5699   31.1102   92.3380   43.0207   18.4816   90.4881   97.9748   43.8870   11.1119   25.8065   40.8720
    %     59.4896   26.2212   60.2843   71.1216   22.1747   11.7418   29.6676   31.8778   42.4167   50.7858    8.5516   26.2482   80.1015    2.9220   92.8854  73.0331   48.8609   57.8525   23.7284   45.8849   96.3089   54.6806   52.1136   23.1594   48.8898   62.4060   67.9136   39.5515   36.7437   98.7982]';
    
    % center = [500 500];
    % rad = 100+400*sqrt(rand(30,1));
    % ang = -pi*2*pi*rand(30,1);
    % xx = rad.*cos(ang)+ center(1);
    % yy = rad.*sin(ang)+ center(2);
    % w_30 = [xx yy]';
    % scatter(w_30(1,:),w_30(2,:))
    
    w_30 = [  806.4430  302.3063  748.1467   22.7443  779.2214  822.2208  447.0217  146.6408  459.5589  415.3604  299.1728  813.4738   96.3286  172.7098  148.1598  191.5577  657.0483  641.6717  946.8172  558.6510  352.9196  408.2349  254.8337  401.3407  348.9908  860.6927  552.3631  679.3186  913.8124   92.0821
        657.7272  112.3110  390.4004  398.9871  418.6168  246.1412  343.4483  215.8258  183.8423  172.6632  657.7292  691.3844  737.8629  148.9883  577.9559  491.9223   93.4039   905.4794  573.0099  846.6488  571.1598  983.1848  783.8461  965.9233  947.7831  803.2735  842.2420  795.0221  496.4085  378.0237]';
    w_30 = w_30/10;
    
    
    
    
    if (op == 100)
        w = w_30(1:N,:);
        beta_o = 10^(-85/10);
    else
        w = 10*w_30(1:N,:);
        beta_o = 10^(-50/10);
    end
    % D_30 = 1e5+1e5*rand(N,1);
    D_30 = 1e4 * [   1.1556    1.2643    1.0828    1.3010    1.1315    1.3270    1.3446    1.3741    1.2253    1.0419    1.1145    1.4567    1.0762    1.4129    1.2692 ...
        1.4981    1.0391    1.2213    1.0533    1.4809    1.0023    1.3875    1.4087    1.4343    1.0422    1.1999    1.1299    1.4000    1.2157    1.4553 ]' ;
    D = D_30(1:N);
    O = 2*D;
    f_max = 5e9; % cycles/sec
    % F_30 = 500+300*rand(1,30);
    F_30 = [651.7871  709.7230  767.2710  787.7874  664.1647  541.5873  544.7882  577.2525  752.2152  576.2847  744.2854  573.0575  778.7791  604.9951  558.9786 ...
        575.3252  684.8134  641.9867  605.4979  749.2486  675.5792  664.9171  775.1581  585.7517  727.1601  726.1187  614.1338  670.3465  522.7563  516.1850]'; % cycles/bit
    F = F_30(1:N);
    c_max = 60e3;
    B = 180e3;
    BW = B;
    
    sf = 1e5;
    
    f_max = f_max/sf;
    D = D/sf;
    O = O/sf;
    c_max = c_max/sf;
    B = B/sf;
    kap = kap*sf^3;
    
    a1 = 2.3; % path-loss exponent
    a2 = 2.5;

    
    
    pc = 0.1; pu = 0.1;
    p_max = 5;
    Gam = 10^(8.2/10); % SNR gap
    K1 = 0.01;
    K2 = 0.99;
    K3 = -4.7;
    K4 = 8.9;
    
    noise_power_desnity = 10.^(-174/10)*1e-3;
    s2 = noise_power_desnity*BW;
    gam = pu*beta_o/(s2*Gam);
    gam_b = pc*beta_o/(s2*Gam);
    
    iter_vec = [];
    
    %% initialization
    % q = zeros(M,2); % UAV horizontal position
    q_ini = [20 20 ; 80 80 ; 50 50];
    if (op == 100)
        q = q_ini(1:M,:);
    else
        q = 10*q_ini(1:M,:);
    end
    f = f_max/N*ones(N,M);
    c = zeros(N,M);
    iter = 0;
    
    %% Association optimiztion
    % d = zeros(N,M);
    
    loop_exit = 0;
    obj_old = 1000;
    while (loop_exit == 0)
        iter = iter + 1;
        
        
        d = pdist2([w zeros(N,1)],[q H]);
        db = pdist2([csX csY 0],[q H])';
        v = zeros(N,M);
        vb = zeros(M,1);
        for j = 1:M
            v(:,j) = H(j)./d(:,j);
            vb(j) = H(j)./db(j);
        end
        % r = zeros(N,M);
        % for i = 1:N
        %     for j = 1:M
        %         r(i,j) = B*log2(1+(K1 + K2/(1+exp(-(K3+K4*v(i,j)))))*gam/d(i,j)^alpha_o);
        %     end
        % end
        r = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*v)))).*gam./(d.^a1));
        rb = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*vb)))).*gam./(db.^a2));
        f = f_max/N*ones(N,M); % because by saying f = f_hat, we will limit the scope of alpha to be that dependant upon alpha
        
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = (1-c(i,j))*D(i)/rb(j);
            end
        end
        
        cvx_begin quiet
        variable mut
        variable alpha_h(N,M)
        % variable f(N,M)
        
        minimize (mut)
        subject to
        % f >= 0;
        mut >= 0;
        for i = 1:N
            sum(alpha_h(i,:)) == 1;
            for j = 1:M
                0 <= alpha_h(i,j) <= 1;
            end
        end
        % for j = 1:M
        %     sum(f(:,j)) <= f_max;
        % end
        for j = 1:M
            mut >= sum(alpha_h(:,j).*(tu(:,j)+tc(:,j) + ts(:,j)));
        end
        cvx_end
        alpha = alpha_h > 0.5;
        
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        mu_vec = zeros(M,1);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = (1-c(i,j))*D(i)/rb(j);
            end
            mu_vec(j) = sum(alpha(:,j).*(tu(:,j)+tc(:,j) + ts(:,j)));
        end
        
        if (dispi == 1)
            display(['iter: ' num2str(iter) ' "First opt" '  ' mu: ' num2str(max(mu_vec)) ' status: ' num2str(cvx_status) ]);
        end
        
        %% Caching policy optimization
        
        
        mu_c = zeros(N,M);
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        cons_c = zeros(M,1);
        c = zeros(N,M);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = D(i)/rb(j);
                            mu_c(i,j) = alpha(i,j).*(tu(i,j)+tc(i,j) + ts(i,j));
%                 mu_c(i,j) = alpha(i,j).*(ts(i,j));
            end
            %         mu_c(mu_c(:,j)==0,j)=inf;
            [sorti, indi] = sort(mu_c(:,j),'descend');
            continue_loop = 1;
            step = 1;
            while (continue_loop == 1)
                c(indi(step),j) = 1;
                cons_c(j) = cons_c(j) + c(indi(step),j)*D(indi(step));
                if (cons_c(j) > c_max)
                    c(indi(step),j) = 0;
                    continue_loop = 0;
                end
                step = step+1;
            end
        end
        
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        mu_vec = zeros(M,1);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = (1-c(i,j))*D(i)/rb(j);
            end
            mu_vec(j) = sum(alpha(:,j).*(tu(:,j)+tc(:,j) + ts(:,j)));
        end
        
        if (dispi == 1)
            display(['iter: ' num2str(iter) ' "Second opt" '  ' mu: ' num2str(max(mu_vec)) ' status: ' num2str(cvx_status) ]);
        end
        
        %% f optimization
        
        
        
        cvx_begin quiet
        variable mut
        variables f(N,M) fb(N,M)
        expressions mu_RHS(N,M) pow_cons(M,1)
        
        for j = 1:M
            pow_cons(j) = pu;
            for i = 1:N
                pow_cons(j) = pow_cons(j) + alpha(i,j)*pow_pos(f(i,j),3)*kap;
                mu_RHS(i,j) = alpha(i,j)*(O(i)/r(i,j) + D(i)*F(i)*fb(i,j)  + (1-c(i,j))*D(i)/rb(j));
            end
        end
        
        minimize (mut)
        subject to
        f >= 0; fb >= 0;
        mut >= 0;
        for j = 1:M
            mut >= sum(mu_RHS(:,j));
            sum(f(:,j).*alpha(:,j)) <= f_max;
            pow_cons(j) <= p_max;
            
            for i = 1:N
                if (alpha(i,j) == 0) % following Y. Zhou paper
                    f(i,j) <= 1e-3;
                end
                [fb(i,j) 1 ; 1 f(i,j)] == semidefinite(2);
            end
        end
        
        cvx_end
        
        
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        mu_vec = zeros(M,1);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = (1-c(i,j))*D(i)/rb(j);
            end
            mu_vec(j) = sum(alpha(:,j).*(tu(:,j)+tc(:,j) + ts(:,j)));
            %             mu_vec(j) = sum(alpha(:,j).*(tu(:,j)+tc(:,j)) );
        end
        if (dispi == 1)
            display(['iter: ' num2str(iter) ' "Third opt" ' ' mu: ' num2str(max(mu_vec)) ' mut: ' num2str(mut) ' ' num2str(cvx_status)]);
        end
        
        
        
        %% Horizontal Position
        
        
        q_hat = q;
        d_hat = pdist2([w zeros(N,1)],[q_hat H]);
        db_hat = pdist2([csX csY 0],[q_hat H])';
        
        cvx_begin quiet
        variable mut
        variables q(M,2) v(N,M) z(N,M) y(N,M) zb(N,M) %f(N,M) fb(N,M)
        expressions r_lb(N,M) v_lb(N,M) mu_RHS(N,M) pow_cons(M,1)
        variables vb(M,1) zz(M,1) yb(M,1) zzb(M,1)
        expressions rb_lb(M,1) vb_lb(M,1)
        
        %         fb = 1./f;
        v_hat = zeros(N,M);
        vb_hat = zeros(M,1);
        for j = 1:M
            v_hat(:,j) = H(j)./d_hat(:,j);
            vb_hat(j) = H(j)/db_hat(j);
        end
        r_hat = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*v_hat)))).*gam./(d_hat.^a1));
        rb_hat = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*vb_hat)))).*gam./(db_hat.^a2));
        for j = 1:M
            pow_cons(j) = pu;
            
            c1 = K1*(1+exp(-(K3+K4*vb_hat(j))))+K2;
            c2 = (1+exp(-(K3+K4*vb_hat(j))))*(db_hat(j)^a2);
            den1 = log(2)*(1+exp(-(K3+K4*vb_hat(j))))*(gam*c1+c2);
            den2 = log(4)*db_hat(j)^2*(gam*c1+c2);
            psiXb = gam*K2*B/den1;
            psiYb = gam*a2*B*c1/den2;
            
            rb_lb(j) = rb_hat(j) - psiXb*(exp(-(K3+K4*vb(j)))-exp(-(K3+K4*vb_hat(j)))) ...
                - psiYb*(yb(j)^2 - norm(q_hat(j,:)-[csX csY])^2);
            den3 = 2*(norm(q_hat(j,:)-[csX csY])^2+H(j)^2)^1.5;
            vb_lb(j) = vb_hat(j) - H(j)/den3*(yb(j)^2 - norm(q_hat(j,:)-[csX csY])^2);
            
            
            for i = 1:N
                pow_cons(j) = pow_cons(j) + alpha(i,j)*pow_pos(f(i,j),3)*kap;
                c1 = K1*(1+exp(-(K3+K4*v_hat(i,j))))+K2;
                c2 = (1+exp(-(K3+K4*v_hat(i,j))))*(d_hat(i,j)^a1);
                den1 = log(2)*(1+exp(-(K3+K4*v_hat(i,j))))*(gam*c1+c2);
                den2 = log(4)*d_hat(i,j)^2*(gam*c1+c2);
                psiX = gam*K2*B/den1;
                psiY = gam*a1*B*c1/den2;
                
                r_lb(i,j) = r_hat(i,j) - psiX*(exp(-(K3+K4*v(i,j)))-exp(-(K3+K4*v_hat(i,j)))) ...
                    - psiY*(y(i,j)^2 - norm(q_hat(j,:)-w(i,:))^2);
                den3 = 2*(norm(q_hat(j,:)-w(i,:))^2+H(j)^2)^1.5;
                v_lb(i,j) = v_hat(i,j) - H(j)/den3*(y(i,j)^2 - norm(q_hat(j,:)-w(i,:))^2);
                
                mu_RHS(i,j) = alpha(i,j)*(O(i)*zb(i,j) + D(i)*F(i)/f(i,j)  + (1-c(i,j))*D(i)*zzb(j));
            end
        end
        
        minimize (mut)
        subject to
        %         f >= 0; fb >= 0;
        mut >= 0;
        y >= 0; z >= 0; v >= 0; zb >= 0;
        yb >= 0; zz >= 0; vb >= 0; zzb >= 0;
        for j = 1:M
            Xmin <= q(j,1); q(j,1) <= Xmax; Ymin <= q(j,2); q(j,2) <= Ymax;
            mut >= sum(mu_RHS(:,j));
            %             sum(f(:,j).*alpha(:,j)) <= f_max;
            %             pow_cons(j) <= p_max;
            
            [zzb(j) 1 ; 1 zz(j)] == semidefinite(2);
            norm(q(j,:)-[csX csY]) <= yb(j);
            vb(j) <= vb_lb(j);
            zz(j) <= rb_lb(j);
            
            
            for i = 1:N
                %                 if (alpha(i,j) == 0) % following Y. Zhou paper
                %                     f(i,j) <= 1e-3;
                %                 end
                [zb(i,j) 1 ; 1 z(i,j)] == semidefinite(2);
                %                 [fb(i,j) 1 ; 1 f(i,j)] == semidefinite(2);
                norm(q(j,:)-w(i,:)) <= y(i,j);
                v(i,j) <= v_lb(i,j);
                z(i,j) <= r_lb(i,j);
            end
        end
        
        cvx_end
        
        d = pdist2([w zeros(N,1)],[q H]);
        db = pdist2([csX csY 0],[q H])';
        v = zeros(N,M);
        vb = zeros(M,1);
        for j = 1:M
            v(:,j) = H(j)./d(:,j);
            vb(j) = H(j)./db(j);
        end
        r = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*v)))).*gam./(d.^a1));
        rb = B*log2(1+(K1 + K2./(1+exp(-(K3+K4*vb)))).*gam./(db.^a2));
        
        tu = zeros(N,M);
        tc = zeros(N,M);
        ts = zeros(N,M);
        mu_vec = zeros(M,1);
        for j = 1:M
            for i = 1:N
                tu(i,j) = O(i)/r(i,j);
                tc(i,j) = D(i)*F(i)/f(i,j);
                ts(i,j) = (1-c(i,j))*D(i)/rb(j);
            end
            mu_vec(j) = sum(alpha(:,j).*(tu(:,j)+tc(:,j) + ts(:,j)));
        end
        if (dispi == 1)
            display(['iter: ' num2str(iter) ' "Fourth opt" ' ' mu: ' num2str(max(mu_vec)) ' mut: ' num2str(mut) ' ' num2str(cvx_status)]);
        end
        obj_fun = max(mu_vec);
        if ( abs((obj_fun/1e-3 - obj_old/1e-3)/(obj_fun/1e-3)) < g_tol)
            loop_exit = 1;
        end
        obj_old = obj_fun;
        
        display(['iter: ' num2str(iter)  ' max_mu: ' num2str(max(mu_vec)) ]);
        
        
    end
    
    display(['N: ' num2str(N) ' max_mu: ' num2str(max(mu_vec)) ' iter: ' num2str(iter)]);
    result_vec(N==N_vec) =  max(mu_vec);
end

plot(N_vec,result_vec*1e3);

result_vec