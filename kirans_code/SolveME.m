function prob = SolveME(H, dissRate, time, vectorised, init)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dims = length(H(1, :));
    nb_sinks = dims - 7;
    nb_modes = 7;
    id = eye(dims);
    % set up init state
    rho0 = zeros(size(H)); 
    rho0(init, init) = 1;


    if vectorised % vectorise rho
        
        % construct lindblad term
        LindbladTerm = sqrt(dissRate)*...
            diag([ones(1, nb_modes), zeros(1, nb_sinks)]);
        
        % write H and L in vectorised notation
        % Note this is transpose even if L is non-hermetian
        preH = kron(H, id); 
        postH = kron(id, transpose(H));
        preDiss = kron(LindbladTerm, id);
        postDiss = kron(id, transpose(LindbladTerm));

        
        base = -1i*(preH - postH) ...
            + preDiss*postDiss - 1/2*(preDiss^2 + postDiss^2);
        S = expm(base*time);
        rhov = reshape(rho0,[dims^2, 1]);
        
        rhot = S*rhov;
        rhot = reshape(rhot, size(rho0));
        entropy = rhot*rhot;
        disp(['tr rho^2 = ', num2str(abs(trace(entropy)))])
        prob = abs(diag(rhot));
        
    else
        if dissRate > 0
            dist('warning matrix version not allow dissRate>0')
        end

        U = expm(-1i*H*time);
        Udag = ctranspose(U);

        Rho = U*rho0*Udag;
        prob = diag(Rho);

    end
end


