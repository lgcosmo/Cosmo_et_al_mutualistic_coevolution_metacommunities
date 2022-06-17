function rbind_3d(;A::AbstractArray)
    #Perform vertical concatenation of the 3rd dimension of an array
    newA = vcat([A[:,:,i] for i in 1:size(A,3)]...)
    return(newA)
end

function rbind_2d(;A::AbstractArray)
    #Perform vertical concatenation of the 3rd dimension of an array
    newA = vcat([A[i,:] for i in 1:size(A,1)]...)
    return(newA)
end

function paste0(;x::String, y::Any)
    r=collect(y)
    result=repeat([x], size(r,1)).*string.(r)
    return result
end


function moore_neighborhood(;n_patches, n, periodic=true) #Create spatial network with moore neighborhood
    
    G=zeros(n_patches, n_patches)
    
    for j in 1:n
        for l in 1:n
            for i in 1:n
                for k in 1:n

                    id_current=((j-1)*n)+i
                    id_nb=((l-1)*n)+k

                    if (abs(k-i) <= 1) && (abs(l-j) <= 1)
                        G[id_current, id_nb]=1.0
                    end

                    if periodic==true
                        if i==1
                            up=((j-1)*n)+n
                            G[id_current, up]=1.0
                        end
                        if j==1
                            left=((n-1)*n)+i
                            G[id_current, left]=1.0
                        end
                        if i==n
                            down=((j-1)*n)+1
                            G[id_current, down]=1.0
                        end
                        if j==n
                            right=((1-1)*n)+i
                            G[id_current, right]=1.0
                        end
                        if j==1 && (i!=1 && i != n)
                            topleft=((n-1)*n)+(i-1)
                            bottomleft=((n-1)*n)+(i+1)
                            G[id_current, topleft]=1.0
                            G[id_current, bottomleft]=1.0
                        end
                        if j==n && (i!=1 && i != n)
                            topright=((1-1)*n)+(i-1)
                            bottomright=((1-1)*n)+(i+1)
                            G[id_current, topright]=1.0
                            G[id_current, bottomright]=1.0
                        end
                        if i==1 && (j!=1 && j != n)
                            tleft=(((j-1)-1)*n)+n
                            tright=(((j+1)-1)*n)+n
                            G[id_current, tleft]=1.0
                            G[id_current, tright]=1.0
                        end
                        if i==n && (j!=1 && j != n)
                            bleft=(((j-1)-1)*n)+1
                            bright=(((j+1)-1)*n)+1
                            G[id_current, bleft]=1.0
                            G[id_current, bright]=1.0
                        end
                        if i==1 && j==1
                            tl11=((n-1)*n)+n
                            tr11=(((j+1)-1)*n)+n
                            lb11=((n-1)*n)+(i+1)
                            G[id_current, tl11]=1.0
                            G[id_current, tr11]=1.0
                            G[id_current, lb11]=1.0
                        end
                        if i==1 && j==n
                            tr1n=((1-1)*n)+n
                            tl1n=(((j-1)-1)*n)+n
                            br1n=((1-1)*n)+(i+1)
                            G[id_current, tr1n]=1.0
                            G[id_current, tl1n]=1.0
                            G[id_current, br1n]=1.0
                        end
                        if i==n && j==1
                            bln1=((n-1)*n)+1
                            tln1=((n-1)*n)+(i-1)
                            brn1=(((j+1)-1)*n)+1
                            G[id_current, bln1]=1.0
                            G[id_current, tln1]=1.0
                            G[id_current, brn1]=1.0
                        end
                        if i==n && j==n
                            brnn=((1-1)*n)+1
                            trnn=((1-1)*n)+(i-1)
                            blnn=(((j-1)-1)*n)+1
                            G[id_current, brnn]=1.0
                            G[id_current, trnn]=1.0
                            G[id_current, blnn]=1.0
                        end
                    end


                end
            end
        end
    end

    G[diagind(G)].=0.0

    return G

end

function occupancy(;z::Array{Float64}) #Compute average occupancy

    O=zeros(size(z,3), size(z,2))

    @inbounds for t in 1:size(z,3)      
       @inbounds for i in 1:size(z,2)

            occ=0.0

           @inbounds for k in 1:size(z,1)

                if z[k,i,t] != 0.0
                    occ+=1.0
                end
            end

            O[t,i]=occ/size(z,1)
            
        end
    end

    return(O)

end

function sp_extinctions(;O::Array{Float64}, min_occup::Float64)
    
    ext=zeros(size(O,2))
    time_ext=zeros(size(O,2))

  @inbounds for i in 1:size(O,2)

        @views O_i=O[:,i]

        if any(O_i.<=min_occup)
            ext[i]+=1.0
            time_ext[i]=findfirst(x-> x<=min_occup, O_i)
        end
    end

    return (ext, time_ext)

end

function local_tm(;z::Array{Float64}, sp_type::Array{String}, α::Float64)

    TM=zeros(size(z,1), size(z,2), size(z,3))
    N_TM=zeros(size(z,1), size(z,2), size(z,3))
    L_TM=zeros(size(z,3), size(z,2))

   @inbounds for t in 1:size(z,3)
    @inbounds for i in 1:size(z,2)
        @inbounds for j in 1:size(z,2)

                if (i==j) || sp_type[i]==sp_type[j] #Skip iteration if j=s, i.e. intraspecific interaction
                    continue
                end

                @inbounds for k in 1:size(z,1)

                    if z[k,i,t]==0.0 || z[k,j,t]==0.0
                        continue
                    end

                    TM[k,i,t]=TM[k,i,t]+exp(-α*(z[k,j,t]-z[k,i,t])^2)
                    N_TM[k,i,t]=N_TM[k,i,t]+1.0
                end
            end
        end
    end

    @. TM=TM/N_TM

    @inbounds for t in 1:size(TM, 3)
        @inbounds for i in 1:size(TM, 2)

            @views M=TM[:,i,t]
            L_TM[t,i]=mean(filter(x-> x>0.0, M))

        end
    end 

    return L_TM

end

function regional_tm(;z::Array{Float64}, sp_type::Array{String}, α::Float64)

    TM=zeros(size(z,1), size(z,2), size(z,3))
    N_TM=zeros(size(z,1), size(z,2), size(z,3))
    R_TM=zeros(size(z,3), size(z,2))

    @inbounds for t in 1:size(z,3)
        @inbounds for i in 1:size(z,2)
            @inbounds for j in 1:size(z,2)

                if (i==j) || sp_type[i]==sp_type[j] #Skip iteration if j=s, i.e. intraspecific interaction
                    continue
                end

                @inbounds for k in 1:size(z,1)

                    if z[k,i,t]==0.0
                        continue
                    end

                    @inbounds for p in 1:size(z,1)

                        if (k==p) || z[p,j,t]==0.0
                            continue
                        end

                        TM[k,i,t]=TM[k,i,t]+exp(-α*(z[p,j,t]-z[k,i,t])^2)
                        N_TM[k,i,t]=N_TM[k,i,t]+1.0

                    end
                end
            end
        end
    end

    @. TM=TM/N_TM

    @inbounds for t in 1:size(TM, 3)
        @inbounds for i in 1:size(TM, 2)

            @views M=TM[:,i,t]
            R_TM[t,i]=mean(filter(x-> x>0.0, M))

        end
    end 

    return R_TM

end

function local_em(;z::Array{Float64}, θ::Array{Float64}, α::Float64)

    EM=zeros(size(z,1), size(z,2), size(z,3))
    L_EM=zeros(size(z,3), size(z,2))

   @inbounds for t in 1:size(z,3)
        @inbounds for i in 1:size(z,2)
            @inbounds for k in 1:size(z,1)

                if z[k,i,t]==0.0
                    continue
                end

                EM[k,i,t]=exp(-α*(θ[k,i,t]-z[k,i,t])^2)

            end
        end
    end

    @inbounds for t in 1:size(EM, 3)
        @inbounds for i in 1:size(EM, 2)

            @views M=EM[:,i,t]
            L_EM[t,i]=mean(filter(x-> x>0.0, M))

        end
    end 

    return L_EM

end

function theta_change(;θ::Array{Float64}, time_ext::Array{Float64}, climchange::Float64)

    C=zeros(size(θ,1), size(θ,2))

    @inbounds for i in 1:size(θ,2)
        @inbounds for k in 1:size(θ,1)

            C[k,i]=((time_ext[i]*climchange)+θ[k,i,1])/θ[k,i,1]
        end
    end

    C_mean=dropdims(mean(C, dims=1), dims=1)

    return C_mean

end









