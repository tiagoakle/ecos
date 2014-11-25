%Runs a large experiment where the iteration count vs 
%the complexity is plotted for different random problems.
% All problems are 
%primal-dual feasible.
%Those with exponential and socp cones are formed of groups of 2 exponential cones and one socp cone
% of size 6. 
%Those of only socp cones are formed of 4 cones of size 3
%We grow m like n and keep it like 1/10;

rand(0);
randn(0);

samples_per_size = 1;
sizes            = [10:10:300];

results = zeros(length(sizes)*samples_per_size,5);
rix = 1;

for size_ix = 1:length(sizes) 
     cones = sizes(size_ix);
     for sample = 1:samples_per_size
      n = 12*cones;
      m = ceil(n/10);
      A = randn(m,n);
      A = sprand(m,n,0.3);
      A = sparse(A);
      G = speye(n);

      %Generate the socp problem exampe 
      s_0 = randn(12*cones,1);
   
      for j=1:3:12*cones
        s_0(j) = norm(s_0(j+1:j+2))+1;
      end
     
      z_0 = randn(12*cones,1);
      for j=1:3:12*cones
        z_0(j) = norm(z_0(j+1:j+2))+1;
      end
       
      x_0 = randn(n,1);
      y_0 = randn(m,1);
      b = A*x_0;
      h = G*x_0 + s_0; 
      c = A'*y_0+G'*z_0;
 
    %Now call xopt
    dims = struct;
    dims.l = 0;
    dims.q = 3*ones(cones*4,1);
    dims.e = 0;
    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.verbose = 1;
    [x,y,info,s,z]=ecos(c,G,h,dims,A,b,opts);
    
    socp_info = info;

    %Prepare the exponential cone problem
    s_0 = randn(12*cones,1); 
    for j=1:6:6*cones
      s_0(j) = norm(s_0(j+1:j+5))+1;
    end
    for j=6*cones+1:3:12*cones
        s_0(j:j+2) = [-1.051383945322714;
                      1.258967884768947;
                      0.556409619469370];
    end
 
    z_0 = randn(12*cones,1);
    for j=1:6:6*cones
      z_0(j) = norm(z_0(j+1:j+5))+1;
    end
    for j=6*cones+1:3:12*cones
      z_0(j:j+2) = [-1.051383945322714;
                      1.258967884768947;
                      0.556409619469370];
    end
     
    x_0 = randn(n,1);
    y_0 = randn(m,1);
    b = A*x_0;
    h = G*x_0 + s_0; 
    c = A'*y_0+G'*z_0;
 
    %Now call xopt
    dims = struct;
    dims.l = 0;
    dims.q = 6*ones(cones,1);
    dims.e = 2*cones;
    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.verbose = 1;
    opts.centrality = 1-log(2); 
    [x,y,info_mixed,s,z]=ecos(c,G,h,dims,A,b,opts);
    
    %Run again with no centrality 
    dims = struct;
    dims.l = 0;
    dims.q = 6*ones(cones,1);
    dims.e = 2*cones;
    opts   = struct;
    opts.maxit = 300;
    opts.feastol = 1.e-6;
    opts.abstol  = 1.e-7;
    opts.verbose = 1;
    opts.centrality = 1e10;
    [x,y,info_mixed_no_cent,s,z]=ecos(c,G,h,dims,A,b,opts);

    %Prepare the exponential cone problem    
    results(rix,:) = [cones*8,cones*12,socp_info.iter,info_mixed.iter,info_mixed_no_cent.iter];
    rix = rix+1;
    end
end

csvwrite('SocpVsMixedUnbounded.csv',results);
%Extract the results and plot

