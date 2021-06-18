function [ZH_vals, ZH_vecs, ZH_inv_vecs] = get_comp_vecs_algo_multi(ws, dt, nD)
    
    dc_vec = gen_dc(nD);
    ZH_vals = [];
    ZH_vecs = [dc_vec];
    for ii = 1:length(ws)
        ZH_vals = [ZH_vals; gen_vals(ws(ii),dt)];
        ZH_vecs = [ZH_vecs gen_vecs(ws(ii),dt,nD)];   
    end  
    ZH_inv_vecs=pinv(ZH_vecs);
    % remove dc vec
    ZH_vecs = ZH_vecs(:,2:end);
    ZH_inv_vecs = ZH_inv_vecs(2:end,:);

end

function ZH_vals= gen_vals(w,dt)  
    cont_vals = [w*i;-w*i];
    ZH_vals = exp(cont_vals*dt);
end

function dc_vec = gen_dc(nD);
    mags = sqrt(1/(nD+1));% normalize magnitude of vector = 1
    dc_vec = ones(nD+1,1)*mags;
end

function ZH_vecs = gen_vecs(w,dt,nD)
    mags = sqrt(1/(nD+1));% normalize magnitude of vector = 1
    del_angle = w*dt;% change in angle each step
    if mod(nD+1,2)==0 % even number of states in delay embed
        start = pi-del_angle/2;
    else
        start = pi;
    end
    angles = [start];
    for jj = 1:ceil((nD+1)/2)-1
        angles = [angles; angles(end)-del_angle];
    end  
    if mod(nD+1,2)==0 % even number of states in delay embed
        angles = [-flipud(angles);angles];
    else
        angles = [-flipud(angles(2:end));angles];
    end        
    vec = mags.*exp(angles*i);% convert mag and angle to complex num
    ZH_vecs = [vec,real(vec)-imag(vec)*i];
end