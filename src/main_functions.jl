#---------------#
#-Colonizations-#
#---------------#

function migrate!(;z::AbstractArray, θ::AbstractArray, G::AbstractArray, sp_type::AbstractArray, mig::AbstractArray, migs::AbstractArray, mig_attempts::AbstractArray, mig_success::AbstractArray, flow::Float64, ρ::Float64, mi::Float64, time::Int)

    #Perform colonizations, modify the array z at time t

    fill!(mig, 0.0) #Resetting migrations
    fill!(migs, 0.0) #Resetting migrations

    @inbounds for s in 1:size(z,2) #Loop through all s species
        @inbounds for i in 1:size(z,1) #Loop through all i sites of species s

            if z[i,s,time]==0.0 #Skip iteration if there is no population of species s at site i
                continue
            end

            @inbounds for j in 1:size(G,2) #Loop through all sites j that can receive colonizations from i
                if G[i,j]==0.0 #If site j cannot receive colonization from i, skip iteration
                    continue
                end

                ksum=0.0 #Sum of potential mutualistic species that habit site j
                tm=0.0 #Variable to store cumulative trait matching with all mutualistic species that habit site j

                @inbounds for k in 1:size(z,2) #Loop through all other k species that may be at site j

                    if (k==s) || (z[j,k,time]==0.0) || sp_type[k]==sp_type[s] #Skip iteration if k=s, i.e. intraspecific interaction, or if there is no population of species k at site j
                        continue
                    end
        
                    tm=tm+(z[j,k,time]-z[i,s,time])^2 #Calculating cumulative trait difference of species s at site i with all other k mutualists at site j
                    ksum=ksum+1.0
                end

                meantm=ifelse(ksum>0.0, tm/ksum, 0.0)
                mut_suit=exp(-ρ*(mi*meantm)) #Mean trait matching of pop i of species s at site j
                env_suit=exp(-ρ*((1.0-mi)*((θ[j,s,(time-1)]-z[i,s,time])^2))) #Environmental matching of pop i of species s at site j in the presence of mutualists
                env_suit2=exp(-ρ*(((θ[j,s,(time-1)]-z[i,s,time])^2))) #Environmental matching of pop i of species s at site j in the absence of mutualists
                migs[i,j,s]=ifelse(meantm == 0.0, env_suit2, mut_suit*env_suit) #Environmental suitability of pop i of species s at site j

                mig_attempts[j,s,time]=mig_attempts[j,s,time]+1.0 #Quantifying attempts of migrations to site j from species s

                if migs[i,j,s] >= rand() #Migration trial

                    mig[i,j,s]=z[i,s,time] #If succesfull, add z value of pop i of species s to migrants to site j
                    mig_success[j,s,time]=mig_success[j,s,time]+1.0 #Updating sucessfull migration

                end
                
                if z[j,s,time]>0.0 && flow>0.0 #If there is already a population of species s at site j
                    mig[i,j,s]=flow*mig[i,j,s] #Population j will only receive a fraction (flow) of population i
                elseif z[j,s,time]>0.0 && flow==0.0
                    mig[i,j,s]=0.0
                end
                
            end
        end
    end

    genflow=dropdims(sum(mig, dims=1), dims=1) #Computing total gene flow to each site (column sums of mig array)
    n_pop=dropdims(sum(x->x>0.0, mig, dims=1), dims=1) #Computing total number of migrant populations to each site

    if flow==0.0 # In the absence of gene flow, only empty sites are colonized and the population with the largest suitability will colonize the site

        for s in 1:size(genflow, 2)
            for j in 1:size(genflow, 1)
               
                if n_pop[j,s]>1.0
                    
                    @views id=findmax(migs[:,j,s])[2]

                    genflow[j,s]=z[id,s,time]
                end
            end
        end
    end

    
    @inbounds for s in 1:size(z,2) #Loop to update trait values after colonization
        @inbounds for i in 1:size(z,1)

            if genflow[i,s]>0.0 #Update trait value if colonization occurred

                z[i,s,time]=ifelse(z[i,s,time]>0.0, (1.0-flow)*z[i,s,time] + genflow[i,s]/n_pop[i,s], genflow[i,s]/n_pop[i,s])

            else

                z[i,s,time]=ifelse(z[i,s,time]>0.0, z[i,s,time], genflow[i,s])
                
            end
        end
    end

end

#-------------#
#-Extinctions-#
#-------------#

function extinct!(;z::AbstractArray, ST::AbstractArray, ext_attempts::AbstractArray, ext_success::AbstractArray, time::Int)
    
    #Perform local extinctions, modify the array z at time t

    
    @inbounds for i in 1:size(z,2) #Loop through each i species
        @inbounds for k in 1:size(z,1) #Loop through each population of species i at site k

            if z[k,i,time]==0.0 #Skip iteration if there is no population of species i at site k
                continue
            end

            ext_attempts[k,i,time]=ext_attempts[k,i,time]+1.0 #Extinction attempt

            if (ST[k,i,time]) <= rand() #Population become extinct if trial is succesfull
                z[k,i,time]=0.0 #Setting population as extinct
                ext_success[k,i,time]=ext_success[k,i,time]+1.0 #If extinction sucesfull, add 1 to extinctions
            end

        end
    end

end

#--------------------------#
#-Networks of interactions-#
#--------------------------#

function interactions!(;z::Array{Float64}, sp_type::Array{String}, A::Array{Float64}, α::Float64, time::Int)
    
    fill!(A, 0.0) #Resetting interactions

    @inbounds for k in 1:size(z,1) #Loop through each k site
        @inbounds for i in 1:size(z,2) #Loop through each i species

            if z[k,i,time]==0.0 #If there is no population of species i at site k, skip iteration
                continue
            end

            @inbounds for j in 1:size(z,2) #Loop through all potential j partners of species i at site k

                if (i==j) || (z[k,j,time]==0.0) || sp_type[i]==sp_type[j] #Avoid intraspecific interactions, interactions of species of the same set, no populations of species j
                    continue
                end

                intprob=(exp(-α*(z[k,j,time]-z[k,i,time])^2)) #Probability of interaction with j and i at site k, proportional to trait matching
                                
                if intprob >= rand() #Interaction trial
                    A[i,j,k]=1.0
                end

                if A[i,j,k]==1.0 #If species i interact with j, force j to interact with i (symmetric interactions)
                    A[j,i,k]=1.0
                end

            end
        end

    end

end

#---------------------------#
#-Evolution and coevolution-#
#---------------------------#

function evolve!(;z::Array{Float64}, A::Array{Float64}, Q::Array{Float64}, θ::Array{Float64}, mi::Float64, α::Float64, ρ::Float64, σ::Float64, time::Int)

    fill!(Q, 0.0) #Resetting Q-matrix

    @inbounds for k in 1:size(z,1) # Loop trough all k populations
        @inbounds for i in 1:size(z,2) # Look through all i species
            
            if z[k,i,time]==0.0 #Skip if there is no population of species i at patch kk
                continue
            end

            @inbounds for j in 1:size(z,2) #Loop through all possible partners of i at site k

                if (i==j) || (z[k,j,time]==0.0) || (A[i,j,k]==0.0) #Skip if there is no population of partner j or i do not interact with j
                    continue
                end

                Q[i,j,k]=exp(-α*(z[k,j,time]-z[k,i,time])^2) #Computing qij

            end
        end
    end

    Q_sum=dropdims(sum(Q, dims=2), dims=2)
    
    @inbounds for k in 1:size(Q,3)
        @inbounds for j in 1:size(Q,2)
            @inbounds for i in 1:size(Q,1)
                Q[i,j,k]=ifelse(Q_sum[i,k]>0.0, Q[i,j,k]/Q_sum[i,k], 0.0)
                Q[i,j,k]=Q[i,j,k]*(z[k,j,time]-z[k,i,time])
            end
        end
    end

    mut=permutedims(dropdims(sum(Q, dims=2), dims=2), (2,1)) #Summing rows of matrix Q and transposing to be at the same order as z array

    @inbounds for i in 1:size(z,2)
        @inbounds for k in 1:size(z,1)

            if z[k,i,time]==0.0
                continue
            end

            z[k,i,(time+1)] = ifelse(mut[k,i]>0.0, z[k,i,time] + σ*ρ*(mi*mut[k,i] + (1.0-mi)*(θ[k,i,time]-z[k,i,time])), z[k,i,time] + σ*ρ*(θ[k,i,time]-z[k,i,time]))
        end
    end
end

#---------------------------------------#
#-Initial populations and traits values-#
#---------------------------------------#

function init_pop!(;z::AbstractArray, time::Int)

    area=trunc(Int, sqrt(size(z,1)))

    @inbounds for i in 1:size(z,2)

        id=sample(1:size(z,1), area , replace=false)

        z[id,i,time].=rand(0.01:0.01:10.0, length(id))

    end
    
end

#--------------------------------------#
#-Initial theta values for all patches-#
#--------------------------------------#

function theta_init!(;θ)
    theta_init=rand(0.01:0.01:10.0, size(θ,1), size(θ,2))
    @inbounds for t in 1:size(θ,3)
        @views θ[:,:,t].=theta_init
    end
end

#----------------#
#-Climate change-#
#----------------#

function clim_change!(;θ::AbstractArray, change::Float64, time::Int64)

    @inbounds for j in 1:size(θ, 2)
        @inbounds for i in 1:size(θ, 1)
            θ[i,j,(time+1)]=θ[i,j,time]+change
        end
    end
end
    

#-------------------------------------#
#-Suitability and Mean trait matching-#
#-------------------------------------#

function suitability!(;z::Array{Float64}, θ::Array{Float64}, ST::Array{Float64}, A::Array{Float64}, ρ::Float64, mi::Float64, time::Int)
    
    #Time index for species environmental optima with environmental change

    if time>1
        time_θ=time-1
    else
        time_θ=time
    end

    @inbounds for k in 1:size(z,1) #Loop through each k site
        @inbounds for i in 1:size(z,2) #Loop through each species i at site k

            if z[k,i,time]==0.0 #Skip iteration if there is no population of species i at site k and set TM and ST as 0
                ST[k,i,time]=0.0
                continue
            end

            tm=0.0 #Quantifying cumulative trait matching
            jsum=0.0 #Quantifying number of interacting partners

            @inbounds for j in 1:size(z,2) #Loop through all j possible partners at site k
                if (i==j) || (z[k,j,time]==0.0) || (A[i,j,k]==0.0) #Skip intraspecific interactions, non-interacting species and absent of populations of j species
                    continue
                end

                tm=(z[k,j,time]-z[k,i,time])^2 #Calculating cumulative trait matching of species i at site k with all other j mutualists at site j
                jsum=jsum+1.0
            end

            meantm=ifelse(jsum>0.0, tm/jsum, 0.0)
            mut_suit=exp(-ρ*(mi*meantm))
            env_suit=exp(-ρ*((1.0-mi)*((θ[k,i,time_θ]-z[k,i,time])^2)))
            env_suit2=exp(-ρ*(((θ[k,i,time_θ]-z[k,i,time])^2))) #Environmental matching of pop k of species i at site k
            ST[k,i,time]=ifelse(meantm == 0.0, env_suit2, mut_suit*env_suit) #Suitability
            
        end
    end
end

#---------------#
#-Main function-#
#---------------#

function coevo_metacom(;n_sp::Int, G::Array{Float64}, climchange::Float64, mi::Float64, α::Float64, ρ::Float64, σ::Float64, flow::Float64, tmax::Int, sim::Int)
    
    #Initializing model
    
    n_p=Int(n_sp/2) #Number of plants
    n_a=Int(n_sp/2) #Number of animals
    sp_type=vcat(repeat(["p"], n_p), repeat(["a"], n_a)) #Vector to classify species as plants or animals
    z=zeros(size(G,1), n_sp, tmax) #Array to store trait values
    ST=zeros(size(G,1), n_sp, tmax) #Array to store suitability values
    A=zeros(size(z,2), size(z,2), size(z,1)) #Array to store adjacency matrices at each time step
    Q=zeros(size(z,2), size(z,2), size(z,1)) #Array to store Q matrices at each time step
    M=zeros(size(G,1), size(G,2), size(z,2)) #Array to store colonization matrices at each time step
    MS=zeros(size(G,1), size(G,2), size(z,2)) #Array to store colonization matrices at each time step
    theta=zeros(size(G,1), n_sp, tmax) #Array to store theta values at each time step
    ext_attempts=zeros(size(G,1), n_sp, tmax) #Array to store extinction attempts at each time step
    ext_success=zeros(size(G,1), n_sp, tmax) #Array to store succesfull extinctions at each time step
    mig_attempts=zeros(size(G,1), n_sp, tmax) #Array to store colonization attempts at each time step
    mig_success=zeros(size(G,1), n_sp, tmax) #Array to store succesfull colonizations at each time step

    theta_init!(θ=theta) #Setting initial theta values
    init_pop!(z=z, time=1) #Setting initial population trait values
    interactions!(z=z, sp_type=sp_type, A=A, α=α, time=1) #Setting initial interaction matrices at each patch
    suitability!(z=z, θ=theta, ST=ST, A=A, ρ=ρ, mi=mi, time=1) #Setting initial suitability

    for t in 1:(tmax-1)

        evolve!(z=z, A=A, Q=Q, θ=theta, α=α, mi=mi, ρ=ρ, σ=σ, time=t) #Evolution and coevolution

        migrate!(z=z, G=G, θ=theta, sp_type=sp_type, mig=M, migs=MS, mig_attempts=mig_attempts, mig_success=mig_success, ρ=ρ, flow=flow, mi=mi, time=t+1) #Colonizations

        interactions!(z=z, sp_type=sp_type, A=A, α=α, time=t+1) #Recalculating interactions after migrations

        suitability!(z=z, θ=theta, ST=ST, A=A, ρ=ρ, mi=mi, time=t+1) #Calculating suitability after migrations

        extinct!(z=z, ST=ST, ext_attempts=ext_attempts, ext_success=ext_success, time=t+1) #Extinctions

        interactions!(z=z, sp_type=sp_type, A=A, α=α, time=t+1) #Recalculating interactions after extinctions

        if climchange>0.0 #Ending simulation if all species underwent extinct

            @views if all(x-> x==0.0, z[:,:,t])
                break
            end
        end

        clim_change!(θ=theta, change=climchange, time=t) #Climate change

    end

    ext_ratio=dropdims(sum(ext_success, dims=1), dims=1)./dropdims(sum(ext_attempts, dims=1), dims=1) #Computing rate of extinction
    ext_ratio=permutedims(ext_ratio, (2,1))
    mig_ratio=dropdims(sum(mig_success, dims=1), dims=1)./dropdims(sum(mig_attempts, dims=1), dims=1) #Computing rate of colonization
    mig_ratio=permutedims(mig_ratio, (2,1))

    occup=occupancy(z=z) #Computing patch occupancy

    sp_e=sp_extinctions(O=occup, min_occup=0.0) #Computing species extinction

    t_change=theta_change(θ=theta, time_ext=sp_e[2], climchange=climchange) #Total change in theta

    ltm=stack(DataFrame(local_tm(z=z, sp_type=sp_type, α=α), :auto)) #Computing local trait matching
    rtm=stack(DataFrame(regional_tm(z=z, sp_type=sp_type, α=α), :auto)) #Computing local regional trait matching
    lem=stack(DataFrame(local_em(z=z, θ=theta, α=α), :auto)) #Computing local environmental matching

    ext_ratio=stack(DataFrame(ext_ratio, :auto))
    mig_ratio=stack(DataFrame(mig_ratio, :auto))

    results=DataFrame(occup, :auto) #Data frame with results
    
    #Formatting dataframe with results
    rename!(results, [(Symbol("x$i")=>Symbol("SP$i")) for i in 1:n_sp])
    results[!, :time]=1:tmax
    results=stack(results)
    results[!,:sp_type]=repeat(sp_type, inner=tmax)
    select!(results, [:time, :variable, :sp_type, :value])
    rename!(results, :variable=>:species,:value=>:occupancy)
    results[!, :local_tm]=ltm.value
    results[!, :local_em]=lem.value
    results[!, :regional_tm]=rtm.value
    results[!, :ext_ratio]=ext_ratio.value
    results[!, :mig_ratio]=mig_ratio.value
    results[!, :sp_extinct]=repeat(sp_e[1], inner=tmax)
    results[!, :time_extinct]=repeat(sp_e[2], inner=tmax)
    results[!, :change_extinct]=repeat(t_change, inner=tmax)

    @. results[!,:n_sp]=n_sp
    @. results[!,:n_a]=n_a
    @. results[!,:n_p]=n_p
    @. results[!,:mi]=mi
    @. results[!,:flow]=flow
    @. results[!,:alpha]=α
    @. results[!,:gvar]=σ
    @. results[!,:rho]=ρ
    @. results[!,:climchange]=climchange
    @. results[!,:simulation]=sim

    return(results)

end