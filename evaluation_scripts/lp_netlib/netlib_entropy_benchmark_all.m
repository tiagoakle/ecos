% Solves the etnropy maximization problem
% using gm scale for scaling and then stores the norm of the solutions to rescale the problem
% min sum_i x_ilogx_i Ax = b

  clear all
  addpath '../'
  %Load the name for the result file
  result_file = getenv('RESULT_FILE');
  if(isempty(result_file))
      result_file ='PDCO_ECOS_Scaled_NScaled_1_20.mat';
  end

  %Load the file that contains the indices for the
  %ufget netlib lps which are in standard form
  load 'standard_form_indices.mat'; 
  problem_count = length(st_ix);
  tuned_scalings = struct;
  tuned_scalings.lp_25fv47.gz = 1; 
  tuned_scalings.lp_25fv47.gx = 1; 

  tuned_scalings.lp_bnl1.gy = 1.e2;
  tuned_scalings.lp_bnl1.gz = 1.e2;

  tuned_scalings.lp_cre_a.gy =1.e2;
  tuned_scalings.lp_cre_a.gz =1.e2;

  %original x 300 y 1000 s 300 z 1000
  tuned_scalings.lp_cre_a.gx =1.e2;
  tuned_scalings.lp_cre_a.gy =1.e2;
  tuned_scalings.lp_cre_a.gz =sqrt(1.e3);

  tuned_scalings.lp_cre_b.gx =1;
  tuned_scalings.lp_cre_b.gy =1.e2;
  tuned_scalings.lp_cre_b.gz =1.e2;
  tuned_scalings.lp_cre_b.potential = true; 

  tuned_scalings.lp_cre_c.gx =100;
  tuned_scalings.lp_cre_c.gy =1000;
  tuned_scalings.lp_cre_c.gz =1000;

  tuned_scalings.lp_cre_d.gx =100;
  tuned_scalings.lp_cre_d.gy =1000;
  tuned_scalings.lp_cre_d.gz =5;
  tuned_scalings.lp_cre_d.gt =1/10;

  %Choose a problem from the list
  for problem_index = 1:problem_count
    
    %Extract the problem 
    problem_uf_ix = st_ix(problem_index);
    %Get the problem from ufget
    P = UFget(problem_uf_ix);
     
    %P = UFget(654);
    %Extract the name
    prob_name = [P.name];
    %Substitute front slash for space
    prob_name(find(prob_name=='/'))=' ';
    
    %Extract the problem data and build the problem structure
    %Problem parameters
    
    [m,n] = size(P.A);    

    %Scale the problem      
    [cscale, rscale] = gmscale([P.A P.b],0.9,100);
    cscale = cscale(1:end-1);
    P.A = diag(sparse(1./rscale))*P.A*diag(sparse(1./cscale));
    P.b = diag(sparse(1./rscale))*P.b;
 
    %Run the experiment using pdco and no scalings
    %---------------------------------------------
 
    pdObj =  get_pdcoEntropy(cscale);
   
    pdOptions = pdcoSet;
    x0    = ones(n,1);
    y0    = zeros(m,1);
    z0    = ones(n,1);
    xsize = 1;                 % Estimate of norm(x,inf) at solution
    zsize = 1;                 % Estimate of norm(z,inf) at solution

    options = pdcoSet;
    options.mu0       = 1e-0;  % An absolute value
    options.LSMRatol1 = 1e-9;  % LSMR and MINRES must solve quite accurately
    options.FeaTol    = 1e-6;
    options.OptTol    = 1e-7;
    options.wait      = 1;     % Allow options to be reviewed before solve
    options.MaxIter   = 300;
    options.wait      = 0;
    
    [pdx,pdy,pdz,pdinform,PDitns,CGitns,time] = ...
    pdco(pdObj,P.A,P.b,zeros(n,1),ones(n,1)*1e10,1.e-6,1.e-6,options,x0,y0,z0,xsize,zsize);
    pd_obj_val = sum((pdx./cscale).*log(pdx./cscale));
    pd_lin_feas = norm(rscale.*(P.A*(pdx./cscale)-P.b));

    problem_r = struct;
    problem_r.('nx') = max(abs(pdx));
    problem_r.('ny') = max(abs(pdy));
    problem_r.('nz') = max(abs(pdz));
    problem_r.('linres') = pd_lin_feas;
    problem_r.('obj')    = pd_obj_val;
    problem_r.('iter')   = PDitns;
    problem_r.('flag')   = pdinform;
    problem_results.(P.name(10:end)).pdco_unscaled = problem_r;

    %Run again using PDCO and setting xize and zsize
    %-----------------------------------------------   
    xsize = problem_r.('nx');
    zsize = problem_r.('nz');
    [pdx,pdy,pdz,pdinform,PDitns,CGitns,time] = ...
    pdco(pdObj,P.A,P.b,zeros(n,1),ones(n,1)*1e10,1.e-6,1.e-6,options,x0,y0,z0,xsize,zsize);
    pd_obj_val = sum((pdx./cscale).*log(pdx./cscale));
    pd_lin_feas = norm(rscale.*(P.A*(pdx./cscale)-P.b));

    problem_r = struct;
    problem_r.('nx') = max(abs(pdx));
    problem_r.('ny') = max(abs(pdy));
    problem_r.('nz') = max(abs(pdz));
    problem_r.('linres') = pd_lin_feas;
    problem_r.('obj')    = pd_obj_val;
    problem_r.('iter')   = PDitns;
    problem_r.('flag')   = pdinform;
    problem_r.('scalings') = struct;
    problem_r.('scalings').xsize = xsize;
    problem_r.('scalings').zsize = zsize;
    problem_results.(P.name(10:end)).pdco_scaled = problem_r;
 

    %Build the conic representation and run ecos with no scalings
    %------------------------------------------------------------

    % min sum(xlog(x)) st Ax=b
    % using u,v,w and wexp(u/w)<v => u < wlog(v/w) => u<-wlog(w/v) => wlog(w/v)<-u
    % minimize -u 0 0  
    % create the constraints Aw=b, v = 1 u,v,w \in K_exp
    % [0 I 0] [u]   = 1
    % [0 0 A] [v]     b
    %         [w]
    
    %Permute!
    %    [-I    ][u]       [su]  
    %[P] [  -I  ][v] + [P] [sv] = 0 
    %    [    -I][w]       [sw] 
    
   
    %Build the conic representation
    ref = zeros(3*n,1);
    ref(1:3:end) = [1:n];
    ref(2:3:end) = [n+1:2*n];
    ref(3:3:end) = [2*n+1:3*n];
    Pr = [[1:3:3*n]';[2:3:3*n]';[3:3:3*n]'];
    Pc = [[1:n]';[n+1:2*n]';[2*n+1:3*n]'];
    G = sparse(Pr,Pc,-ones(3*n,1));
    h  = zeros(3*n,1);
    cscaled  = [-1./cscale;zeros(n,1);1./cscale];
    c = cscaled;
    bb = [ones(n,1);P.b];
    AA =[ [sparse(n,n) speye(n,n)  sparse(n,n)];
          [sparse(m,n) sparse(m,n)        P.A ]];

    
    %Now call ecos
    dims = struct;
    dims.l = 0;
    dims.q = [];
    dims.e = n;
    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.reltol  = 1.e-7;
    opts.centrality = 1e10;
    opts.potential = false;
    [x,y,info,s,z]=ecos(c,G,h,dims,AA,bb,opts);

   
    problem_r = struct;
    problem_r.('nx') = max(abs(x));
    problem_r.('ny') = max(abs(y));
    problem_r.('ns') = max(abs(s));
    problem_r.('nz') = max(abs(z));
    problem_r.('linres') = norm(rscale.*(P.A*(s(3:3:end)./cscale)-P.b));
    problem_r.('obj')    = sum((s(3:3:end)./cscale).*log(s(3:3:end)./cscale));
    problem_r.('iter')   = info.iter;
    problem_r.('flag')   = info.exitflag;
    problem_results.(P.name(10:end)).ecos_unscaled = problem_r;

    %Repeat the solution using ECOS with scalings 
    %---------------------------------------------------------
    %Scale the conic representation
    nx = max(abs(x));
    ns = max(abs(s));
    nz = max(abs(z));
    ny = max(abs(y));
    
    gx = sqrt(max(nx,1));
    gz = sqrt(max(nz,1))/sqrt(max(ns,1)); 
    gy = sqrt(max(ny,1));
    gt = 1;
  
    %Check if there are hand-tuned scalings for any of 
    %these problems
    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.reltol  = 1.e-7;
    opts.centrality = 1e10;
    opts.potential = false;

    tuned_scaling_flag = false;
    if(isfield(tuned_scalings,P.name(10:end)))
        scalings = tuned_scalings.(P.name(10:end));
        if(isfield(scalings,'gx')) gx = scalings.gx; end
        if(isfield(scalings,'gy')) gy = scalings.gy; end
        if(isfield(scalings,'gz')) gz = scalings.gz; end
        if(isfield(scalings,'gt')) gt = scalings.gt; end
        if(isfield(scalings,'potential')) opts.potential = scalings.potential; end
        tuned_scaling_flag = true;
    end
    %Now use the scaling on the problem
    AA = AA*gx*gy;
    G  = G*gz*gx;
    bb  = bb*gy*gt; 
    c  = c*gx*gt;
    h  = h*gz*gt;
    
    %Now call xopt
    dims = struct;
    dims.l = 0;
    dims.q = [];
    dims.e = n;
    [x,y,info,s,z]=ecos(c,G,h,dims,AA,bb,opts);

    problem_r = struct;
    problem_r.('nx') = max(abs(x));
    problem_r.('ny') = max(abs(y));
    problem_r.('ns') = max(abs(s));
    problem_r.('nz') = max(abs(z));
    problem_r.('linres') = norm(rscale.*(P.A*(s(3:3:end)./cscale)-P.b));
    problem_r.('obj')    = sum((s(3:3:end)./cscale).*log(s(3:3:end)./cscale));
    problem_r.('iter')   = info.iter;
    problem_r.('flag')   = info.exitflag;
    problem_r.('scalings') = struct;
    problem_r.('scalings').('gx') = gx;
    problem_r.('scalings').('gz') = gz;
    problem_r.('scalings').('gy') = gy;
    problem_r.('scalings').('gt') = gt; 
    problem_r.('scalings').('tuned') = tuned_scaling_flag;
    problem_results.(P.name(10:end)).ecos_scaled = problem_r;
    

    %Now call the MOSEK solver ------------------- 
    param = struct;
    param.MSK_DPAR_INTPNT_NL_TOL_PFEAS = 1.e-6;
    param.MSK_DPAR_INTPNT_NL_TOL_DFEAS = 1.e-7;
    param.MSK_DPAR_INTPNT_NL_TOL_REL_GAP = 1.e-7;
    d = 1./cscale;
    c = -1./cscale.*(log(cscale));
 
    diary on
    [mskres] = mskenopt(d,c,P.A,P.b,P.b,param)
    diary off
    %Load the diary and extract the number of iterations..... what a hack :(
    fid = fopen('diary');
    diary_string = fread(fid);
    pattern = 'Interior-point          - iterations : (\d*)';
    itrcount = regexp(char(diary_string'),pattern,'tokens');
    fclose(fid);
    %Delete the diary 
    delete('./diary')

    %Calculate the norm of the residual, and the objective value
    problem_r = struct;
    sol_nz        = find(mskres.sol.itr.xx > 0);
    problem_r.('linres') = norm(rscale.*(P.A*(mskres.sol.itr.xx./cscale)-P.b));
    problem_r.('obj')    = sum(mskres.sol.itr.xx(sol_nz)./cscale(sol_nz)).*log(mskres.sol.itr.xx(sol_nz)./cscale(sol_nz));
    problem_r.('iter')   = itrcount{1}{1}; 
    problem_r.('flag')   = mskres.rcode;
    problem_results.(P.name(10:end)).mosek_scaled = problem_r;

    %For good measure save the results up until now
    tuned_scalings_used = tuned_scalings;
    save(result_file,'problem_results','tuned_scalings_used');
end

field_names = fieldnames(problem_results);
%Print the results 
fprintf( 'problem nm:      PU  FG| PS  FG| EU  FG| ES  FG| MS  FG\n');
for j = 1:length(field_names)
    field_name = field_names(j);
    problem_result = problem_results.(field_name{:});
    pui = problem_result.pdco_unscaled.iter;
    puf = problem_result.pdco_unscaled.flag;
    psi = problem_result.pdco_scaled.iter;
    psf = problem_result.pdco_scaled.flag;
    eui = problem_result.ecos_unscaled.iter;
    euf = problem_result.ecos_unscaled.flag;
    esi = problem_result.ecos_scaled.iter;
    esf = problem_result.ecos_scaled.flag;
    msi = problem_result.mosek_scaled.iter;
    msf = problem_result.mosek_scaled.flag;
  
    fprintf('%15s %3i %2i |%3i %2i |%3i %2i |%3i %2i |%3s %2i |\n',field_name{:},pui,puf,psi,psf,eui,euf,esi,esf,msi,msf);
end
