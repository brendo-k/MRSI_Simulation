% Modified code from Ilya Kuprov from spiniach package!
%
% Calculates exponential propagator exp(-1i*L*timestep) using scaled 
% and squared Taylor series method. Syntax:
%
%                 P=propagator(spin_system,L,timestep)
%
% Parameters:
%
%    L          -  Hamiltonian or Liouvillian matrix
%
%    timestep   -  propagation time step
%
% Outputs:
%
%    P          -  propagator matrix
%
% Note: GPUs are supported, add 'gpu' to sys.enable array during 
%       calculation setup.
%
% Note: propagator caching (https://doi.org/10.1063/1.4928978) is
%       supported, add 'caching' to sys.enable array during calcu-
%       lation setup.
%
% Note: we did have Chebyshev and Newton series here at one point,
%       as well as the Pade method. None of them lived up to their
%       marketing.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk
%
% <http://spindynamics.org/wiki/index.php?title=propagator.m>

function P=propagator(spin_system,L,timestep)

% Check consistency
grumble(L,timestep);

% Set a shorthand for -i*L*dt
A=-1i*L*timestep;

% Fast bypass for small matrices
if size(A,1) < 200
    P=expm(full(A)); return;
end
% 
% % Check the cache
% if ismember('caching',spin_system.sys.enable)
%     
%     % Generate the cache record name
%     filename=[spin_system.sys.scratch '/p_' md5_hash(A) '.mat'];
%     
%     % Try loading the cache record
%     if exist(filename,'file')
%         report(spin_system,'cache record found and used.'); load(filename,'P'); 
%         report(spin_system,['propagator dimension ' num2str(size(P,1)) ', nnz ' num2str(nnz(P))...
%                             ', density ' num2str(100*nnz(P)/numel(P)) '%, sparsity ' num2str(issparse(P))]);
%         return;
%     else
%         report(spin_system,'cache record not found, computing...'); tic;
%     end
%     
% end

% Estimate the norm
mat_norm=cheap_norm(A);

% Check the norm
if mat_norm>1e9
    
    % The user is doing something silly, bomb out
    error('norm of -i*L*dt exceeds 1e9, check your L and your dt.');

elseif mat_norm>1024
    
    % If the user is really pushing it, take precautionary measures
    report(spin_system,'WARNING - the time step requested greatly exceeds the timescale of system dynamics.');
    report(spin_system,'WARNING - exponentiation tolerance will be set to 1e-14.');
    spin_system.tols.prop_chop=1e-14;
    
elseif mat_norm>16
    
    % Inform the user just in case
    report(spin_system,'WARNING - the time step requested exceeds the timescale of system dynamics.');
    
end

% Determine scaling and squaring parameters
n_squarings=max([0 ceil(log2(mat_norm))]); 
scaling_factor=2^n_squarings;
report(spin_system,['scaling -i*L*dt down by ' num2str(scaling_factor) ...
                    ' and squaring the propagator ' num2str(n_squarings) ' times.']);

% Scale and clean up the matrix
if scaling_factor>1, A=(1/scaling_factor)*A; end
A=clean_up(spin_system,A,spin_system.tols.prop_chop);

% Get the propagator
if ismember('gpu',spin_system.sys.enable)&&(size(A,1)>500)
    
    % Run Taylor series procedure on the GPU
    A=gpuArray(A); 
    P=speye(size(A));
    next_term=gpuArray.speye(size(A));
    n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+gather(next_term); n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on GPU in ' num2str(n) ' iterations.']);
    
else
    
    % Run Taylor series procedure on the CPU
    P=speye(size(A)); next_term=P; n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+next_term; n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on CPU in ' num2str(n) ' iterations.']);
    
end

% Reclaim memory
clear('A','next_term');

% Clean up the result
P=clean_up(spin_system,P,spin_system.tols.prop_chop);

% Inform the user
report(spin_system,['propagator dimension ' num2str(size(P,1)) ', nnz ' num2str(nnz(P))...
                    ', density ' num2str(100*nnz(P)/numel(P)) '%, sparsity ' num2str(issparse(P))]);

% Run the squaring stage
if n_squarings>0
    
    % Inform the user
    report(spin_system,'setting up the squaring process...');
    
    % Run the appropriate squaring process
    if ismember('gpu',spin_system.sys.enable)&&(size(P,1)>500)
        
        % Move the array to GPU
        P=gpuArray(P);
        
        % Run the squaring on the GPU
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['GPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
        
        % Gather the propagator
        P=gather(P);
        
    elseif (~isworkernode)&&(nnz(P)>1e6)&&issparse(P)
        
        % Distribute the propagator
        P=distributed(P);
        
        % Run codistributed CPU squaring
        spmd
            
            % Run the squaring stage
            for n=1:n_squarings
                
                % Inform the user
                if labindex==numlabs
                    report(spin_system,['codistributed CPU squaring step ' num2str(n) '...']);
                end
                
                % Square the propagator
                P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
                
            end
            
        end
        
        % Gather the propagator
        P=gather(P);
        
    else
        
        % Run serial CPU squaring
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['CPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
      
    end
    
    % Inform the user
    report(spin_system,['propagator dimension ' num2str(size(P,1)) ', nnz ' num2str(nnz(P))...
                        ', density ' num2str(100*nnz(P)/numel(P)) '%, sparsity ' num2str(issparse(P))]);

end
    
% Write the cache record
if ismember('caching',spin_system.sys.enable)&&(toc>1)
    
    % Save the propagator
    save(filename,'P','-v7.3'); report(spin_system,'cache record saved.');
    
elseif ismember('caching',spin_system.sys.enable)
    
    % Inform the user
    report(spin_system,'cache record not worth saving.');
    
end

end


% To preserve one's mind intact through a modern college education is a
% test of courage and endurance, but the battle is worth it and the sta-
% kes are the highest possible to man: the survival of reason.
%
% Ayn Rand

