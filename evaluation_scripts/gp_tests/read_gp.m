function [c,G,h,A,b,dims,num_ter,num_con] = read_gp(file_name) 
    %Reads a GP in the format provided by Erling Andersen's .eo files.
    %Sets up the constraints in the form solved by coneopt:
    % 
    %
    % ecos_A = []
    % G = [  , Ix,  ][x]  [si]   [1     ]
    %     [-A,   ,  ][v] +[su] = [log(c)]
    %     [  ,-I ,  ]     [sv]   [0     ]
    %     [  ,   ,  ]     [sw]   [1     ]
    %
    % num_ter: is the total number of monomial terms in the cosnstraints and objective.
    % Each "term" gives rise to a triplet u,v,w which is constrained to be in an exponential cone
    % num_var: is the number of variables in the problem
    % num_con: is the number of posinomial constraints in the problem
    %
    % The matrix G is constructed as follows
    % A is of size num_ter by num_var
    % IX is of size num_con by num_ter, each row corresponds to 
    % a constraint and each column to a term, the corresponding entry 
    % is 1 if the term belongs to the constraint.
    % ix_0 contains 1s in the indices that correspond to the terms in the objective.
   
    %Load the data from the file
    f         = fopen(file_name,'r');
    %Read the number of constraints
    num_con   = fscanf(f,'%i',1);
    %Read the number of variables
    num_var   = fscanf(f,'%i',1);
    %Read number of terms 
    num_ter   = fscanf(f,'%i',1);
    %Read the 'c' coefficients for the terms
    [c_coef,c_coef_count] = fscanf(f,'%g',num_ter);
    if(c_coef_count ~= num_ter); error('Unable to read c coefficients'); end
    
    %Read the constraint indices
    [constraint_index,constraint_index_count] = fscanf(f,'%i',num_ter);
    if(constraint_index_count~=num_ter); error('Unable to read the constraint indices'); end
    constraint_index = constraint_index+1; %The data is zero indexed

    %Concatenate into IJV format, not an efficient implementation
    I = [];
    J = [];
    V = []; 
    
    %Define the vectors that indicate which terms belong to each constraint
    Ix = sparse(constraint_index,[1:num_ter]',ones(num_ter,1),num_con+1,num_ter,num_ter);
    Ix0 = Ix(1,:);
    Ix  = Ix(2:end,:);

    %Define a vector for the objective
    a_0 = zeros(3*num_ter+num_var,1);
    while true
        row = fscanf(f,'%i %i %e',3);
        if isempty(row)
            break;
        end
        %The data is 0 indexed
        I = [I;row(1)+1];
        J = [J;row(2)+1];
        V = [V;row(3)];
    end
    %All done with the file reading
    fclose(f);

    A   = sparse(I,J,V,num_ter,num_var);
    G = [[sparse(num_con,num_var),Ix                        ];...
         [-A                     ,sparse(num_ter,num_ter)   ];...
         [sparse(num_ter,num_var),-speye(num_ter)           ];
         [sparse(num_ter,num_var),sparse(num_ter,num_ter)   ]];

    %Shuffle to match the ecos ordering
    I = [1:3*num_ter+num_con];
    J = [1:num_con]';
    for(j=1:num_ter)
       J = [J;[num_con+j;num_con+j+num_ter;num_con+j+2*num_ter]];
    end
    
    shuffle = sparse(I,J,ones(3*num_ter+num_con,1));
    %Do the shuffle
    G = shuffle*G;   
    A = [];
    b = [];
    dims = struct;
    dims.l = num_con;
    dims.e = num_ter;

    h  = [ones(num_con,1);log(c_coef);zeros(num_ter,1);ones(num_ter,1)];
    h  = shuffle*h; 
    c = [zeros(num_var,1);full(Ix0)'];
    
end

