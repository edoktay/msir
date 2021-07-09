%TSIR_TEST_SS

%Run tests for matrices from SuiteSparse
%Note: requires ssget MATLAB interface

matids = [907, 1641, 293, 906, 462, 464, 1199, 2338, 253];

fprintf('Running SDQ tests for SuiteSparse matrices\n');
for i = 1:numel(matids)
    Problem = ssget(matids(i));
    A = full(Problem.A);
    b = ones(size(A,2),1);
    
    uf = 1; u =2; ur = 4;
    
    namearr = split(Problem.name,'/');
    name = namearr{2};
    fprintf('Testing matrix %s\n', name);
    
    snbase = strcat('figs/',name,'_');
    
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
end

fprintf('Running HSD tests for SuiteSparse matrices\n');
for i = 1:numel(matids)
    Problem = ssget(matids(i));
    A = full(Problem.A);
    b = ones(size(A,2),1);
    
    uf = 0; u =1; ur = 2;
    
    namearr = split(Problem.name,'/');
    name = namearr{2};
    fprintf('Testing matrix %s\n', name);
    
    snbase = strcat('figs/',name,'_');
    
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
end

fprintf('Running HDQ tests for SuiteSparse matrices\n');
for i = 1:numel(matids)
    Problem = ssget(matids(i));
    A = full(Problem.A);
    b = ones(size(A,2),1);
    
    uf = 0; u = 2; ur = 4;
    
    namearr = split(Problem.name,'/');
    name = namearr{2};
    fprintf('Testing matrix %s\n', name);
    
    snbase = strcat('figs/',name,'_');
    
    
    sir(A,b,uf,u,ur,maxit,strcat(snbase,'SIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    sgmresir(A,b,uf,u,ur,maxit,strcat(snbase,'SGMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
  
    gmresir(A,b,uf,u,ur,maxit,strcat(snbase,'GMRESIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow

    
    tsir(A,b,uf,u,ur,maxit, .5,strcat(snbase,'TSIR_',num2str(uf),num2str(u),num2str(ur),'_1'));
    drawnow
    
end


