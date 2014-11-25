clear all
addpath ../
addpath ./gp
%Runs all GP problems in coneopt
%Saves the results to the table stored in 
% the tex file 'benchmark_results_gp.mat'

%Problems in the paper
problem_names = {...
'beck751.eo',...
'beck752.eo',...
'beck753.eo',...
'bss2.eo',...
'car.eo',...
'demb761.eo',...
'demb762.eo',...
'demb763.eo',...
'demb781.eo',...
'fang88.eo',...
'fiac81a.eo',...
'fiac81b.eo',...
'gptest.eo',...
'jha88.eo',...
'mra01.eo',...
'mra02.eo',...
'rijc781.eo',...
'rijc782.eo',...
'rijc783.eo',...
'rijc784.eo',...
'rijc785.eo',...
'rijc786.eo',...
'rijc787.eo'};

results = {{'name','iter','exit flag'}};
%Now iterate over all problems
problem_count = size(problem_names,2);
for(j = 1:problem_count)
    problem_file_name = problem_names{j};
    problem_file_name = ['./gp/',problem_file_name]
    [c,G,h,A,b,dims,num_ter,num_con] = read_gp(problem_file_name);
    dims.q = [];
    dims.ep = dims.e;
    pfeas = [-1.051383945322714;
              1.258967884768947;
              0.556409619469370];

    %Make the initial points
    x0 = zeros(size(G,2),1);
    y0 = zeros(size(A,1),1);
    s0 = [ones(num_con,1);repmat(pfeas,num_ter,1)];
    z0 = [ones(num_con,1);repmat(pfeas,num_ter,1)];
    tau0 = 1;
    kappa0 = 1;

    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.reltol  = 1.e-7;
    opts.centrality = 1e10;
    opts.potential  = true;

    [x,y,info,s,z]=ecos(c,G,h,dims,opts);
    results = {results{:},{problem_names{j},info.iter,info.exitflag}}; 

end

%Print the first result
res = results{1};
fprintf('%s %s %s \n',res{1},res{2},res{3});
for(j=2:length(results))
    res = results{j};
    fprintf('%20s, %3i, %3i, \n',res{1},res{2},res{3});
end

%%Define the cell array for the results and include the header
%results = {{'Prob name','co','status','lstep','status','lstepn','status',...
%                                      'lstep fo','status','predcor','status'}};
%                                     
%fid = 1;
%problem_count = size(problem_names,2);
%for(j =1:problem_count)
%    
%    fprintf('Will solve problem %s \n',problem_names{j});
%    problem_file_name = problem_names{j};
%    %Add the path to the file
%    problem_file_name = ['./gp/',problem_file_name];
%    [AA,bb,cc,num_ter,num_var,num_con] = read_gp(problem_file_name);
%    
%    %----------------------------------------------------
%    % coneopt call 
%    %---------------------------------------------------
%    % build cone:
%    K.npos = 2*num_var+1;
%    K.npow = 0;
%    K.nexp = num_ter;
%    K.nlog = 0;
%    K      = getbarrpar(K);
%    
%    %Set the parameters
%    pars.n = 3*num_ter + 2*num_var + 1;
%    pars.m = 2*num_ter + num_con+1;
%    pars.echo = 4;
%
%    pars.secord = 1;
%    pars.cnbfgsstps = 3;
%    pars.theta = 0.7;
%    pars.eta   = 0.5;
%    pars.beta  = 0.2;
%   
%    pars.rhoP  = 1e-5;
%    pars.rhoD  = 1e-5;
%    pars.rhoA  = 1e-5;
%    pars.rhoG  = 1e-5;
%    pars.rhoI  = 1e-7;
%    pars.rhoM  = 1e-7;
%
%    
%    % starting point:
%    t00 = 1;
%    up0 = ones(num_var,1);
%    um0 = ones(num_var,1);
%    w00 = -ones(num_ter,1);
%    v00 = ones(num_ter,1);  
%    y00 = 0.5*ones(num_ter,1);
%    
%    v0.x  = [t00;up0;um0;w00;v00;y00];
%    
%    % call to coneopt:
%    R = coneopt(AA,bb,cc,v0,K,pars);
%    cone_kkt = R.dat.nkktsolves;
%    cone_sta = R.status;
%   
%    %--------------------------------------------------------------------------
%    % Solve with nscs
%    %--------------------------------------------------------------------------
%    %Set up call to nscs
%    % starting point:
%    t00 = 1;
%    up0 = ones(num_var,1);
%    um0 = ones(num_var,1);
%    w00 = -ones(num_ter,1);
%    v00 = ones(num_ter,1);  
%    y00 = 0.5*ones(num_ter,1);
%    
%    %Extract the problem data and build the problem structure
%    problem = struct;
%    problem.A = AA;
%    problem.b = bb;
%    problem.c = cc;
%    
%    %Problem parameters
%    problem.m =2*num_ter + num_con+1;
%    problem.n = 3*num_ter + 2*num_var + 1;
%    problem.n_free = 0;
%    problem.n_constrained = 2*num_var+1+3*num_ter;
%    problem.n_pos       = 2*num_var+1;
%    problem.soc_cones   = 0;
%    problem.n_soc_cones = 0;
%    problem.n_sdp_cones = 0;
%    problem.sdp_cones     = 0;
%    problem.n_exp_cones   = num_ter;
%    problem.n_power_cones = 0;
%   
%    %Algorithm parameters
%    pars.solve_second_order = true;
%
%    x0c  = [t00;up0;um0;w00;v00;y00];
%    x0f  = [];
%
%    set_default_pars_nscs;
%    pars.max_iter = 200;
%    nscs 
%    nscs_kkt = state.kkt_solves; 
%    nscs_sta = state.exit_reason;
%     
%    set_default_pars_nscs_long_step;
%    pars.max_iter = 200;
%
%    %--------------------------------------------------------------------------
%    % Solve with nscs long step and no nt scaling
%    %--------------------------------------------------------------------------
%    pars.use_nesterov_todd_scaling = false;
%    nscs_long_step
%    nscs_ls_kkt = state.kkt_solves; 
%    nscs_ls_sta = state.exit_reason;
%      
%    %--------------------------------------------------------------------------
%    % Solve with nscs long step and using nt scaling
%    %-------------------------------------------------------------------------- 
%    pars.use_nesterov_todd_scaling = true;
%    nscs_long_step
%    nscs_lsnt_kkt = state.kkt_solves; 
%    nscs_lsnt_sta = state.exit_reason;
% 
%    %--------------------------------------------------------------------------
%    % Solve with nscs long step and using nt scaling and no second order
%    %-------------------------------------------------------------------------- 
%    pars.use_nesterov_todd_scaling = true;
%    pars.solve_second_order        = false;
%    nscs_long_step
%    nscs_lsntfo_kkt = state.kkt_solves; 
%    nscs_lsntfo_sta = state.exit_reason;
%
%    %--------------------------------------------------------------------------
%    % Save the results 
%    %--------------------------------------------------------------------------
%    problem_result =   {problem_names{j},cone_kkt,cone_sta,...
%                        nscs_ls_kkt,nscs_ls_sta,...
%                        nscs_lsnt_kkt,nscs_lsnt_sta,...
%                        nscs_lsntfo_kkt,nscs_lsntfo_sta,...
%                        nscs_kkt,nscs_sta};
%    results        =   {results{:},problem_result};
%
%    %Print out after every iteration 
%       fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
%                                                  results{1}{1},results{1}{2},...
%                                                  results{1}{3},results{1}{4},...
%                                                  results{1}{5},results{1}{6},...
%                                                  results{1}{7},results{1}{8},...
%                                                  results{1}{9},results{1}{10},...
%                                                  results{1}{11});
%    for j=2:size(results,2)
%       fprintf(fid,'%10s, %3i, %15s, %3i, %15s, %3i, %15s, %3i, %15s, %3i, %15s \n',...
%                                                  results{j}{1},results{j}{2},...
%                                                  results{j}{3},results{j}{4},...
%                                                  results{j}{5},results{j}{6},...
%                                                  results{j}{7},results{j}{8},...
%                                                  results{j}{9},results{j}{10},...
%                                                  results{j}{11});
%
%    end
% 
%end
%
%%Clear all but results 
%clear -REGEXP '^(?!.*?results).*'
%%Save to file
%save 'benchmark_results_gp'
%
%
%
