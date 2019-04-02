function V_U_dp_signs=get_V_U_dp(VU_measures,min_V,min_U) % VU_measures - column vector

    % Check small sizes!!!

    N_frames=length(VU_measures);
    D=[repmat(VU_measures,1,min_V) repmat(1-VU_measures,1,min_U)];
    P=zeros(N_frames,min_V+min_U);
    for N=2:N_frames
        for M=1:min_V+min_U
            inds=Get_valid_indeces(min_V,min_U,M);
            [V,I]=max(D(N-1,inds));
            D(N,M)=D(N,M)+V;
            P(N,M)=inds(I);
        end
    end
    
    V_U_dp_pos=zeros(N_frames,1);
    [V,I]=max(D(end,:));
    V_U_dp_pos(end)=I;
    for N=N_frames:-1:2
        V_U_dp_pos(N-1)=P(N,V_U_dp_pos(N));
    end
    V_U_dp_signs=zeros(size(V_U_dp_pos));
    V_U_dp_signs(V_U_dp_pos<=min_V)=1;
end

function inds=Get_valid_indeces(min_V,min_U,cur_ind)
    if(cur_ind==1)
        inds=min_V+min_U;
        return;
    end
    
    if(cur_ind==min_V)
        inds=[cur_ind-1 cur_ind];
        return;
    end
    
    if(cur_ind==min_V+min_U)
        inds=[cur_ind-1 cur_ind];
        return;
    end
    
    inds=cur_ind-1;    
end