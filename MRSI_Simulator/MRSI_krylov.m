% ADAPTED FROM SPINACH BY ILYA KUPROV
%
% Krylov propagation function. Avoids matrix exponentiation, but can be
% slow. Should be used when the Liouvillian exponential does not fit in-
% to the system memory, but the Liouvillian itself does. Syntax:
%
%     answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)
%
% Arguments:
%
%   L      - the Liouvillian to be used during evolution
%
%   rho    - the initial state vector or a horizontal stack thereof
%
%   output - a string giving the type of evolution that is required
%
%             'final' - returns the final state vector or a horizontal
%                       stack thereof.
%
%             'trajectory' - returns the stack of state vectors giving
%                            the trajectory of the system starting from
%                            rho with the user-specified number of steps
%                            and step length.
%
%             'total'   - returns the integral of the observable trace
%                         from the simulation start to infinity. This
%                         option requires the presence of relaxation.
%
%             'multichannel' - returns the time dynamics of several
%                              observables as rows of a matrix. Note
%                              that destination state screening may be
%                              less efficient when there are multiple
%                              destinations to screen against.
%
%   coil   - the detection state, used when 'observable' is specified as
%            the output option. If 'multichannel' is selected, the coil
%            should contain multiple columns corresponding to individual
%            observable vectors.
%
% Outputs:
%
%      answer - a vector or a matrix, depending on the options set during
%               the call.
%
% Note: this function does not support Hilbert space formalisms.
%
% Note: we initially had a faithful implementation of the Krylov process
%       here - subspace, orthogonalisation, projection, etc., but in all
%       our testing it was much inferior to the reordered Taylor process
%       that is currently implemented below.
%
% i.kuprov@soton.ac.uk
%
% <http://spindynamics.org/wiki/index.php?title=krylov.m>

function answer=MRSI_krylov(spin_system,L,coil,rho,timestep,nsteps,output)

    % Check consistency
    grumble(spin_system,L,coil,rho,timestep,nsteps,output);

    % Scale the state vector
    scaling=max([1 norm(rho, 1)]);
    rho=rho/scaling;

    % Compute the generator matrix
    A=-1i*L*timestep;

    % Upload data to GPU or optimise layout
    if ismember('gpu',spin_system.sys.enable)
        A=gpuArray(A);
        rho=gpuArray(full(rho));
        coil=gpuArray(coil);
        location='GPU';
    else
        location='CPU'; rho=full(rho);
    end

    % Estimate the norm of the generator
    norm_mat=cheap_norm(A);

    % Inform the user about matrix densities
    report(spin_system,['-i*L*dt dim ' num2str(size(A,1)) ...
        ', nnz '       num2str(nnz(A)) ...
        ', norm '      num2str(norm_mat) ...
        ', sparsity '  num2str(issparse(A)) ...
        ', polyadic '  num2str(isa(A,'polyadic'))]);

    % Determine the number of substeps and scale the generator
    nsubsteps=ceil(norm_mat/2); A=(1/nsubsteps)*A;

    % Inform the user about the schedule
    report(spin_system,['taking ' num2str(nsteps) ' Krylov steps with '...
        num2str(nsubsteps) ' substeps each.']);
    if nsubsteps>100
        report(spin_system,'WARNING - too many substeps, consider using evolution() here.');
    end

    % Start feedback timer
    feedback=tic();

    % Decide the output type
    switch output

        case 'final'

            % Loop over steps
            % Loop over substeps
            for k=1:nsubsteps
                % Taylor series
                next_term=1; 
                m=1;
                %Multiply first term by rho
                next_term = next_term * rho;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m);
                    rho=rho+next_term; m=m+1;
                end
                if m>32, warning('loss of accuracy in krylov(), use evolution() instead'); end

            end

            % Assign the answer
            answer=gather(rho);

        case 'trajectory'

            % Preallocate the answer and set the starting point
            answer=zeros([size(rho,1) (nsteps+1)],'like',1i);
            answer(:,1)=gather(rho);

            % Loop over steps
            for n=1:nsteps

                % Loop over substeps
                for k=1:nsubsteps

                    % Taylor series
                    next_term=rho; m=1;
                    while nnz(abs(next_term)>eps)>0
                        next_term=(A*next_term)*(1/m);
                        rho=rho+next_term; m=m+1;
                    end
                    if m>32, warning('loss of accuracy in krylov(), use evolution() instead'); end

                end

                % Assign the answer
                answer(:,n+1)=gather(rho);

                % Inform the user
                if (n==nsteps)||(toc(feedback)>1)
                    report(spin_system,[location ' Krylov step ' num2str(n)...
                        ' out of ' num2str(nsteps) ' done.']);
                    feedback=tic();
                end

            end

        case 'multichannel'

            % Preallocate the answer
            answer=zeros([size(coil,2) (nsteps+1)],'like',1i);

            % Set the initial point
            answer(:,1)=answer(:,1)+gather(coil'*rho);

            % Loop over steps
            for n=1:nsteps

                % Loop over substeps
                for k=1:nsubsteps

                    % Taylor series
                    next_term=rho; m=1;
                    while nnz(abs(next_term)>eps)>0
                        next_term=(A*next_term)*(1/m);
                        rho=rho+next_term; m=m+1;
                    end
                    if m>32, warning('loss of accuracy in krylov(), use evolution() instead'); end

                end

                % Assign the answer
                answer(:,n+1)=gather(coil'*rho);

                % Inform the user
                if (n==nsteps)||(toc(feedback)>1)
                    report(spin_system,[location ' Krylov step ' num2str(n)...
                        ' out of ' num2str(nsteps) ' done.']);
                    feedback=tic();
                end

            end

        otherwise

            error('unknown output option.');

    end

    % Scale the anwer back
    answer=scaling*answer;

end

% Consistency enforcement
function grumble(spin_system,L,coil,rho,timestep,nsteps,output)
    if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
        error('this function only works in Lioville space.');
    end
    if ~isnumeric(L)
        error('Liouvillian must be numeric.');
    end
    if ~isnumeric(coil)
        error('coil argument must be numeric.');
    end
    if (~isnumeric(rho))&&(~iscell(rho))
        error('rho argument must either be numeric or a cell array');
    end
    if ~isnumeric(timestep)
        error('timestep argument must be numeric.');
    end
    if ~isnumeric(nsteps)
        error('nsteps must be numeric.');
    end
    if (~ischar(output))||(~ismember(output,{'observable','final',...
            'trajectory','total','multichannel','refocus'}))
        error('observable argument must be a valid character string.');
    end
end

% Freedom can be achieved and retained only by sober men who
% take humanity as it is, not as humanity should be.
%
% Russell Kirk

