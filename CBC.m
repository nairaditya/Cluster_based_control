function [J,bi] = CBC(K,t,area,force_norm,power_norm,gamma,cl_centroid)
%% Nair et al., Cluster-based feedback control of turbulent post-stall separated flows, JFM 2019.

%% Inputs

% K          : number of clusters
% time       : time vector of data collected
% area       : area of the control surface
% force_norm : normalizing force
% power_norm : normalizing power
% gamma      : tradeoff
% cl_centroid: cluter centroids

%% Outputs

% J          : optimized aerodynamic power
% bi         : optimized simplex

%% Requirements

% Navier-Stokes solver (run_Navier_Stokes)
% outputs at each iteration (cluster_amplitudes and control case forces)
% All outputs must be stored in results

%% Initial control cases

Nb               = K+1;                 % number of control laws for simplx
J               = zeros(Nb,1);          % Total power
bi              = lhsdesign(Nb,K);      % Latin hypercube design
for niter = 1:Nb
    try
        cluster_amp = load(['results/output',num2str(niter),...
            '/control.txt']);
        force       = load(['results/output',num2str(niter),...
            '/force.txt']);
        CD_ctrl     = force(:,2)/force_norm;
        CL_ctrl = force(:,3)/force_norm;
        JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
        JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
        Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
        Jb          = (1/t(end))*trapz(t, ...
            abs(cluster_amp).^3*area)/power_norm;
    catch ME
        bi_temp = bi(niter,:)';
        run_Navier_Stokes(bi_temp,cl_centroid);
    end
    J(niter)    = Ja + Jb;
end

%% Simplex optimization

[J,index] = sort (J);                  % sorting cost function
bi        = bi(index,:);               % sorting corresponding control laws
rho       = 1;                          % rho > 0          (reflection)
xi        = 2;                          % xi  > max(rho, 1)(expansion)
gam       = 0.5;                        % 0 < gam < 1      (contraction)
sig       = 0.5;                        % 0 < sig < 1      (shrinkage)
tolerance = 1.0E-06;                    % tolerance for convergence
converged = 0;                          % check convergence
diverged  = 0;                          % check divergence
maxiter   = 250;                        % maximum number of iterations
n_dim     = Nb-1;                       % total dimension of simplex
count     = 1;

while (~converged && ~diverged)
     %% Reflection
    
    bi_bar = sum(bi(1:n_dim,:))/n_dim;                % midpoint of simplex
    bi_r   = (1+rho)*bi_bar - rho*bi(n_dim+1,:);      % reflection
    bi_temp = bi_r';
    niter = niter + 1;
    try
        cluster_amp = load(['results/output',num2str(niter),...
            '/control.txt']);
        force       = load(['results/output',num2str(niter),...
            '/force.txt']);
        CD_ctrl     = force(:,2)/force_norm;
        CL_ctrl = force(:,3)/force_norm;
        JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
        JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
        Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
        Jb          = (1/t(end))*trapz(t, ...
            abs(cluster_amp).^3*area)/power_norm;
    catch ME
        run_Navier_Stokes(bi_temp,cl_centroid);
        exit;
    end
    J_r = Ja + Jb;                                   % evaluate reflection
    indexr = niter;
    
    if ( J(1) <= J_r && J_r <= J(n_dim))             % accept reflection
        bi(n_dim+1,:) = bi_r;
        J(n_dim+1  ) = J_r;
        index(n_dim+1) = indexr;
                
    elseif ( J_r < J(1) )        
        %% Expansion
        
        bi_e = (1+rho*xi)*bi_bar ...
            - rho*xi*bi(n_dim+1,:);
        bi_temp = bi_e';
        niter = niter+1;      
        try
            cluster_amp = load(['results/output',num2str(niter),...
                '/control.txt']);
            force       = load(['results/output',num2str(niter),...
                '/force.txt']);
            CD_ctrl     = force(:,2)/force_norm;
            CL_ctrl = force(:,3)/force_norm;
            JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
            JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
            Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
            Jb          = (1/t(end))*trapz(t, ...
                abs(cluster_amp).^3*area)/power_norm;
        catch ME
            run_Navier_Stokes(bi_temp,cl_centroid);
            exit;
        end
        J_e = Ja + Jb;              % evaluate expansion
        indexe = niter;
        
        if ( J_e < J_r )            % accept expansion
            bi(n_dim+1,:) = bi_e;
            J(n_dim+1  )  = J_e;
            index(n_dim+1) = indexe;
        else
            bi(n_dim+1,:) = bi_r;
            J(n_dim+1  )  = J_r;
            index(n_dim+1) = indexr;
        end     
        
    elseif ( J(n_dim) <= J_r && J_r < J(n_dim+1))        
        %% Outside contraction
        
        bi_c = (1+rho*gam)*bi_bar - rho*gam*bi(n_dim+1,:);
        bi_temp = bi_c';
        niter = niter+1;      
        try
            cluster_amp = load(['results/output',num2str(niter),...
                '/control.txt']);
            force       = load(['results/output',num2str(niter),...
                '/force.txt']);
            CD_ctrl     = force(:,2)/force_norm;
            CL_ctrl = force(:,3)/force_norm;
            JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
            JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
            Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
            Jb          = (1/t(end))*trapz(t, ...
                abs(cluster_amp).^3*area)/power_norm;
        catch ME
            run_Navier_Stokes(bi_temp,cl_centroid);
            exit;
        end
        J_c = Ja + Jb;              % evaluate contraction
        indexc = niter;
        
        if (J_c <= J_r)             % accept contraction
            bi(n_dim+1,:) = bi_c;
            J(n_dim+1  )  = J_c;
            index(n_dim+1) = indexc;
        else           
            %% Shrink landscape
            
            x = bi;
            [~,n_dim_temp] = size ( x );
            i = 1;
            cluster_amp = load(['results/output',num2str(index(i)),...
                '/control.txt']);
            force       = load(['results/output',num2str(index(i)),...
                '/force.txt']);
            CD_ctrl     = force(:,2)/force_norm;
            CL_ctrl = force(:,3)/force_norm;
            JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
            JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
            Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
            Jb          = (1/t(end))*trapz(t, ...
                abs(cluster_amp).^3*area)/power_norm;           
            f = zeros(n_dim_temp+1,1);
            f(i) = Ja + Jb;
            for i = 2 : n_dim_temp + 1
                x(i,:) = sig * x(i,:) + ( 1.0 - sig ) * x(1,:);
                bi_temp = x(i,:)';
                niter = niter+1;         
                try
                    cluster_amp = load(['results/output',num2str(niter),...
                        '/control.txt']);
                    force       = load(['results/output',num2str(niter),...
                        '/force.txt']);
                    CD_ctrl     = force(:,2)/force_norm;
                    CL_ctrl = force(:,3)/force_norm;
                    JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
                    JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
                    Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
                    Jb          = (1/t(end))*trapz(t, ...
                        abs(cluster_amp).^3*area)/power_norm;
                catch ME
                    run_Navier_Stokes(bi_temp,cl_centroid);
                    exit;
                end
                f(i) = Ja + Jb;
                index(i)= niter;
            end
            bi = x;J = f;
        end             
    else                            
        %% Inside contraction
        
        bi_c = (1-gam)*bi_bar + gam*bi(n_dim+1,:);
        bi_temp = bi_c';
        niter = niter+1;      
        try
            cluster_amp = load(['results/output',num2str(niter),...
                '/control.txt']);
            force       = load(['results/output',num2str(niter),...
                '/force.txt']);
            CD_ctrl     = force(:,2)/force_norm;
            CL_ctrl = force(:,3)/force_norm;
            JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
            JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
            Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
            Jb          = (1/t(end))*trapz(t, ...
                abs(cluster_amp).^3*area)/power_norm;
        catch ME
            run_Navier_Stokes(bi_temp,cl_centroid);
            exit;
        end
        J_c = Ja + Jb;               % evaluate contraction
        indexc = niter;
        
        if (J_c < J(n_dim+1))        % accept contraction
            bi(n_dim+1,:) = bi_c;
            J(n_dim+1  ) = J_c;
            index(n_dim+1) = indexc;
        else            
            %% Shrink landscape
            
            x = bi;
            [~,n_dim_temp] = size ( x );
            i = 1;
            cluster_amp = load(['results/output',num2str(index(i)),...
                '/control.txt']);
            force       = load(['results/output',num2str(index(i)),...
                '/force.txt']);
            CD_ctrl     = force(:,2)/force_norm;
            CL_ctrl = force(:,3)/force_norm;
            JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
            JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
            Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
            Jb          = (1/t(end))*trapz(t, ...
                abs(cluster_amp).^3*area)/power_norm;
            f = zeros(n_dim_temp+1,1);
            f(i) = Ja + Jb;
            for i = 2 : n_dim_temp + 1
                x(i,:) = sig * x(i,:) + ( 1.0 - sig ) * x(1,:);
                bi_temp = x(i,:)';
                niter = niter+1;            
                try
                    cluster_amp = load(['results/output',num2str(niter),...
                        '/control.txt']);
                    force       = load(['results/output',num2str(niter),...
                        '/force.txt']);
                    CD_ctrl     = force(:,2)/force_norm;
                    CL_ctrl = force(:,3)/force_norm;
                    JCD_raw     = (1/t(end))*trapz(t, CD_ctrl);
                    JCL_raw     = (1/t(end))*trapz(t, CL_ctrl);
                    Ja          = gamma*(JCD_raw/(JCL_raw^(3/2)));
                    Jb          = (1/t(end))*trapz(t, ...
                        abs(cluster_amp).^3*area)/power_norm;
                catch ME
                    run_Navier_Stokes(bi_temp,cl_centroid);
                    exit;
                end
                f(i) = Ja + Jb;
                index(i)= niter;
            end
            bi = x;J = f;
        end      
    end
    [J,indx] = sort(J);
    bi = bi(indx,:);
    index = index(indx);
    converged = J(n_dim+1)-J(1) < tolerance;
    diverged  = (maxiter<niter);
    count = count + 1;
end