%TSIR_TEST  

%Runs tests for random dense matrices

n = 100;
maxit = 10;

condnums = [1e1,1e2,1e4,1e5,1e7,1e9,1e11,1e14];


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mode 2 matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run sdq tests
fprintf('Running SDQ tests for mode 2 random matrices\n');
for i = 1:numel(condnums)
    
    fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
    
    rng(1);
    A = gallery('randsvd',n,condnums(i),2);
    b = randn(n,1);
    

    uf = 1; u =2; ur = 4;
    
    snbase = strcat('figs/mode2_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
    
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end


%Run hsd tests
fprintf('Running HSD tests for mode 2 random matrices\n');
for i = 1:numel(condnums)
    
     fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
    
    rng(1);
    A = gallery('randsvd',n,condnums(i),2);
    b = randn(n,1);
    
    
    uf = 0; u =1; ur = 2;
    
    snbase = strcat('figs/mode2_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
 
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
 
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    

    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end

%Run hdq tests
fprintf('Running HDQ tests for mode 2 random matrices\n');
for i = 1:numel(condnums)
    
     fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
    
    rng(1);
    A = gallery('randsvd',n,condnums(i),2);
    b = randn(n,1);
    
    
    uf = 0; u =2; ur = 4;
    
    snbase = strcat('figs/mode2_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
     
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mode 3 matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run sqd tests
fprintf('Running SDQ tests for mode 3 random matrices\n');
for i = 1:numel(condnums)
    
     fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
    
    rng(1);
    A = gallery('randsvd',n,condnums(i),3);
    b = randn(n,1);
    
    
    uf = 1; u =2; ur = 4;
    
    snbase = strcat('figs/mode3_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
     
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
   
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
   
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end


%Run hsd tests
fprintf('Running HSD tests for mode 3 random matrices\n');
for i = 1:numel(condnums)
    
     fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
     
    rng(1);
    A = gallery('randsvd',n,condnums(i),3);
    b = randn(n,1);
    
    
    uf = 0; u =1; ur = 2;
    
    snbase = strcat('figs/mode3_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
 
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end

%Run hdq tests
fprintf('Running HDQ tests for mode 3 random matrices\n');
for i = 1:numel(condnums)
    
     fprintf('\nCondition number 1e%s\n',num2str(log10(condnums(i))));
    
    rng(1);
    A = gallery('randsvd',n,condnums(i),3);
    b = randn(n,1);
    
    
    uf = 0; u =2; ur = 4;
    
    snbase = strcat('figs/mode3_rand_size_100_cond_e',num2str(log10(condnums(i))),'_');
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
 
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
       
end