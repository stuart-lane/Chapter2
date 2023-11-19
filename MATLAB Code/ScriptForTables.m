%%% STUART LANE
%%% 07/06/2023 
%%% PRODUCE TABLES 6.1, 7.1 and 7.2 FROM CHAPTER - EMPIRICAL APPLICATION ON
%%% ELASTICITY OF INTERTEMPORAL SUBSTITUTION, USING DATA FROM YOGO (2004).

clear; 

%%% SAVE RESULTS AT END? SET == 1 IF YES, == 0 IF NO:
SaveResults = 1;

%%

tic

rng(234435,'combRecursive');
     
print = zeros(1,8);
Table6_1a = zeros(1,6);
Table6_1b = zeros(1,6);
Table7_1a = zeros(1,7);
Table7_1b = zeros(1,7);
Table7_2a = zeros(1,7);
Table7_2b = zeros(1,7);

num_sim = 10;

for country_number = 1:11

    % SELECT COUNTRY
    if country_number == 1
       dataset = importdata('Data\AULQ.txt');
    elseif country_number ==2
       dataset = importdata('Data\CANQ.txt');
    elseif country_number ==3
       dataset = importdata('Data\FRQ.txt');
    elseif country_number ==4
       dataset = importdata('Data\GERQ.txt');
    elseif country_number ==5
       dataset = importdata('Data\ITAQ.txt');
    elseif country_number ==6
       dataset = importdata('Data\JAPQ.txt');
    elseif country_number ==7
       dataset = importdata('Data\NTHQ.txt');
    elseif country_number ==8
       dataset = importdata('Data\SWDQ.txt');
    elseif country_number ==9
       dataset = importdata('Data\SWTQ.txt');
    elseif country_number ==10
       dataset = importdata('Data\UKQ.txt');
    else
        dataset = importdata('Data\USAQ.txt');
    end
    
    % SELECT NUMBER OF LAGS DEPENDING ON COUNTRY
    if country_number <=10
        L = 4;
    else
        L = 6;
    end
    
    dataset = dataset.data;

    dat = dataset(:,1);
    dc = dataset(:,6);
    r = dataset(:,7);
    ir = dataset(:,8);
    Z = dataset(:,9:12);

    obs = length(dataset);
    onevec = ones(obs,1);

    M = eye(obs) - onevec*onevec'/obs;
    dc = M*dc;
    r = M*r;
    ir = M*ir;
    Z = M*Z;
    Z = fun_orthonormalise(Z,obs);
    izz = inv(Z'*Z);
    
    %%%%%%%%%%%%%%
    % TABLE 6.1a %
    %%%%%%%%%%%%%%
    
    pi_tsls_ir = Z\dc;
    tsls_ir = (Z*pi_tsls_ir)\ir; 
    u_tsls_ir = ir - dc*tsls_ir;      
    v_tsls_ir = dc - Z*pi_tsls_ir;
    
    rho_tsls_ir = corr(u_tsls_ir,v_tsls_ir);
    cov_tsls_ir = cov(u_tsls_ir,v_tsls_ir);
    cov_tsls_ir = cov_tsls_ir(1,2);
    
    [liml_ir,pi_liml_ir,u_liml_ir,v_liml_ir] = fun_liml(ir,dc,Z,izz,obs,1);
    rho_liml_ir = corr(u_liml_ir,v_liml_ir);
    
    % COMPUTE EFFECTIVE F STATISTIC
    Wir = NeweyWest(v_tsls_ir,Z,L,0);
    trWir = trace(Wir);
    effFir = dc'*Z*Z'*dc/trace(Wir);
    effDOFir = (trWir^2)*(1+2*10)/(trace(Wir'*Wir) + 2*10*trWir*max(eig(Wir)));
    cvir = fun_critval(effDOFir,10);

    J_tsls_ir = fun_kpcal_nw(Z,Z(:,2:4),pi_tsls_ir,u_tsls_ir,L);
    KP_liml_ir = fun_kpcal_nw(Z,Z(:,2:4),pi_liml_ir,u_liml_ir,L);

    %%%%%%%%%%%%%%
    % TABLE 6.1b %
    %%%%%%%%%%%%%%
    
    pi_tsls_dc = Z\ir;
    tsls_dc = (Z*pi_tsls_dc)\dc;
    u_tsls_dc = dc - ir*tsls_dc;
    v_tsls_dc = ir - Z*pi_tsls_dc;
    
    rho_tsls_dc = corr(u_tsls_dc,v_tsls_dc);
    cov_tsls_dc = cov(u_tsls_dc,v_tsls_dc);
    cov_tsls_dc = cov_tsls_dc(1,2);

    [liml_dc,pi_liml_dc,u_liml_dc,v_liml_dc] = fun_liml(dc,ir,Z,izz,obs,1);
    rho_liml_dc = corr(u_liml_dc,v_liml_dc);

    J_tsls_dc = fun_kpcal_nw(Z,Z(:,2:4),pi_tsls_dc,u_tsls_dc,L);
    KP_liml_dc = fun_kpcal_nw(Z,Z(:,2:4),pi_liml_dc,u_liml_dc,L);
    
    % COMPUTE EFFECTIVE F STATISTIC
    Wdc = NeweyWest(v_tsls_dc,Z,L,0);
    trWdc = trace(Wdc);
    effFdc = ir'*Z*Z'*ir/trace(Wdc);
    effDOFdc = (trWdc^2)*(1+2*10)/(trace(Wdc'*Wdc) + 2*10*trWdc*max(eig(Wdc)));
    cvdc = fun_critval(effDOFdc,10);
    
    % PUT VALUES INTO TABLES
      
    Table6_1a(end+1,:) = [effFdc cvdc tsls_dc liml_dc  J_tsls_dc KP_liml_dc]
    Table6_1b(end+1,:) = [effFir cvir tsls_ir liml_ir J_tsls_ir KP_liml_ir]
 
    % SIMULATIONS CALIBRATED TO ABOVE ESTIMATES
    
    sb7_1a = zeros(num_sim,6);
    sb7_1b = zeros(num_sim,6);
    sb7_2a = zeros(num_sim,6);
    sb7_2b = zeros(num_sim,6);
    
    for si = 1:num_sim
        
        %%%%%%%%%%%%%%
        % TABLE 7.1a %
        %%%%%%%%%%%%%%
        
        eps = mvnrnd([0 0],[var(u_tsls_dc) cov_tsls_dc; cov_tsls_dc var(v_tsls_dc)],obs);
        
        z = mvnrnd(mean(Z),cov(Z),obs);        
        z = fun_orthonormalise(z,obs);
        
        u = sqrt(abs(z(:,1))).*eps(:,1);
        v = sqrt(abs(z(:,2))).*eps(:,2);
        
        x = z*pi_tsls_dc + v;
        y = x*tsls_dc + u;
        
        M = eye(obs) - ones(obs,1)*ones(obs,1)'/obs;
        y = M*y;
        x = M*x;
        z = M*z;

        izz = inv(z'*z);
        b2sls = (z*(z\x))\y;
        u2sls = y - x*b2sls;
        pi2sls = z\x;
        v2sls = x - z*(z\x);
        x2sls =  z*pi2sls;

        [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,obs,1);
        
        rhoxu = corr(x,u);

        z2 = z(:,2:end);
        xhat = z*piliml;
        mxz = z2 - xhat*(xhat\z2);
       
        kpsample = fun_kpcal_nw(z,z2,piliml,uliml,0);

        score2sls = fun_kpcal_nw(z,z2,pi2sls,u2sls,0); 
                       
        W = NeweyWest(v2sls,z,0,0);
        trW = trace(W);
        effF = x'*z*z'*x/(trace(W));
        
        sb7_1a(si,:) = [effF rhoxu b2sls bliml score2sls kpsample];
        
        %%%%%%%%%%%%%%
        % TABLE 7.2b &
        %%%%%%%%%%%%%%
        
        b2sls = (z*(z\y))\x;
        u2sls = x - y*b2sls;
        pi2sls = z\y;
        v2sls = y - z*(z\y);
        y2sls =  z*pi2sls;

        [bliml,piliml,uliml,vliml] = fun_liml(x,y,z,izz,obs,1);
        
        rhoyu = corr(y,u);

        z2 = z(:,2:end);
        yhat = z*piliml;
        myz = z2 - yhat*(yhat\z2);
       
        kpsample = fun_kpcal_nw(z,z2,piliml,uliml,0);

        score2sls = fun_kpcal_nw(z,z2,pi2sls,u2sls,0); 
       
        W = NeweyWest(v2sls,z,0,0);
        trW = trace(W);
        effF = y'*z*z'*y/(trace(W));
        
        sb7_2b(si,:) = [effF rhoyu b2sls bliml score2sls kpsample];
        
        %%%%%%%%%%%%%%
        % TABLE 7.1b %
        %%%%%%%%%%%%%%
        
        eps = mvnrnd([0 0],[var(u_tsls_ir) cov_tsls_ir; cov_tsls_ir var(v_tsls_ir)],obs);
        
        z = mvnrnd(mean(Z),cov(Z),obs);        
        z = fun_orthonormalise(z,obs);
        
        u = sqrt(abs(z(:,1))).*eps(:,1);
        v = sqrt(abs(z(:,2))).*eps(:,2);
        
        x = z*pi_tsls_ir + v;
        y = x*tsls_ir + u;
        
        M = eye(obs) - ones(obs,1)*ones(obs,1)'/obs;
        y = M*y;
        x = M*x;
        z = M*z;

        izz = inv(z'*z);
        b2sls = (z*(z\x))\y;
        u2sls = y - x*b2sls;
        pi2sls = z\x;
        v2sls = x - z*(z\x);
        x2sls =  z*pi2sls;

        [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,obs,1);
        
        rhoxu = corr(x,u);

        z2 = z(:,2:end);
        xhat = z*piliml;
        mxz = z2 - xhat*(xhat\z2);
       
        kpsample = fun_kpcal_nw(z,z2,piliml,uliml,0);

        score2sls = fun_kpcal_nw(z,z2,pi2sls,u2sls,0); 
                       
        W = NeweyWest(v2sls,z,0,0);
        trW = trace(W);
        effF = x'*z*z'*x/(trace(W));
        
        sb7_1b(si,:) = [effF rhoxu b2sls bliml score2sls kpsample];
        
        %%%%%%%%%%%%%%
        % TABLE 7.2a %
        %%%%%%%%%%%%%%
        
        b2sls = (z*(z\y))\x;
        u2sls = x - y*b2sls;
        pi2sls = z\y;
        v2sls = y - z*(z\y);
        y2sls =  z*pi2sls;

        [bliml,piliml,uliml,vliml] = fun_liml(x,y,z,izz,obs,1);
        
        rhoyu = corr(y,u);

        z2 = z(:,2:end);
        yhat = z*piliml;
        myz = z2 - yhat*(yhat\z2);
       
        kpsample = fun_kpcal_nw(z,z2,piliml,uliml,0);

        score2sls = fun_kpcal_nw(z,z2,pi2sls,u2sls,0); 
       
        W = NeweyWest(v2sls,z,0,0);
        trW = trace(W);
        effF = y'*z*z'*y/(trace(W));
        
        sb7_2a(si,:) = [effF rhoyu b2sls bliml score2sls kpsample];
        
    end
        
    Table7_1a(end+1,:) = [obs mean(sb7_1a(:,1)) mean(sb7_1a(:,2)) median(sb7_1a(:,3:4)) mean(chi2cdf(sb7_1a(:,5:6),3,'upper')<0.05)];
    Table7_1b(end+1,:) = [obs mean(sb7_1b(:,1)) mean(sb7_1b(:,2)) median(sb7_1b(:,3:4)) mean(chi2cdf(sb7_1b(:,5:6),3,'upper')<0.05)];
    Table7_2a(end+1,:) = [obs mean(sb7_2a(:,1)) mean(sb7_2a(:,2)) median(sb7_2a(:,3:4)) mean(chi2cdf(sb7_2a(:,5:6),3,'upper')<0.05)];
    Table7_2b(end+1,:) = [obs mean(sb7_2b(:,1)) mean(sb7_2b(:,2)) median(sb7_2b(:,3:4)) mean(chi2cdf(sb7_2b(:,5:6),3,'upper')<0.05)];
    
    while country_number < 11
        country_number = country_number + 1; 
    end
    
end

Table6_1a_print = round(Table6_1a(2:end,:),3)
Table6_1b_print = round(Table6_1b(2:end,:),3)
Table7_1a_print = round(Table7_1a(2:end,:),3)
Table7_1b_print = round(Table7_1b(2:end,:),3)
Table7_2a_print = round(Table7_2a(2:end,:),3)
Table7_2b_print = round(Table7_2b(2:end,:),3)

round(Table6_1a(2:end,3),3)

if SaveResults == 1
    mkdir TablesofResults
    writematrix(Table6_1a_print,'TablesofResults\Table6_1aValues.txt');
    writematrix(Table6_1b_print,'TablesofResults\Table6_1bValues.txt');
    writematrix(Table7_1a_print,'TablesofResults\Table7_1aValues.txt');
    writematrix(Table7_1b_print,'TablesofResults\Table7_1bValues.txt');
    writematrix(Table7_2a_print,'TablesofResults\Table7_2aValues.txt');
    writematrix(Table7_2b_print,'TablesofResults\Table7_2bValues.txt');
else
end
mean(Table7_2b_print(:,end-1))
toc

[prctile(Table7_2b_print(:,3),97.5) prctile(Table7_2b_print(:,3),2.5)] 