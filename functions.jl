using LinearAlgebra, Symbolics, SymPy, Latexify, LaTeXStrings, JuMP, Catalyst, Combinatorics, Random, Gurobi

"""
extract_useful_from_crn(rn::ReactionSystem)

Extract useful matrices such as stoichiometric matrix, substoichiometric matrix, complexes stoichiometric matrix, incidence matrix from
the reaction network. Also get number of species, reactions and complexes
"""

function extract_useful_from_crn(rn)
    #check if type is the right one, otherwise raise error
    if !(rn isa Catalyst.ReactionSystem)
        error("Expected a Catalyst.ReactionSystem. Got: $(typeof(rn))")
    end
    S = netstoichmat(rn)#compute stoichiometric matrix
    Sub = substoichmat(rn) #compute substoichiometric matrix
    Y = complexstoichmat(rn) #compute complexes stoichiometric matrix
    I = incidencemat(rn) #compute incidence matrix
    nb_complexes = size(I,1) # get number of complexes
    nb_reactions = size(I,2) # get number of reactions
    nb_species = size(S,1) # get number of species
    return S,Sub,Y,I,nb_complexes,nb_reactions,nb_species
end

"""
find_source_complexes(rn::ReactionSystem)

Order the complexes of a reaction network and return a vector of size nb_complexes where
i^th element corresponds to the reactions have i as source complex.
"""
function find_source_complexes(rn)
    if !(rn isa ReactionSystem)
        error("Input must be a Catalyst.ReactionSystem. Got: $(typeof(rn))")
    end
    I = incidencemat(rn) #compute incidence matrix
    I = ordering_complexes(rn)
    # raise potential errors
    if isempty(I)
        error("The incidence matrix is empty. The reaction network may not contain any reactions.")
    end
    if any(x -> x ∉ (-1:1), I)
        error("Unexpected values in the incidence matrix. Expected only -1, 0, or 1.")
    end
    # loop on the incidence matrix to find columns (reactions) where complex i is a source (reactant)
    source_complexes = []
    for i in 1:size(I, 1)
        indices = findall(x -> x == -1, I[i, :])
        push!(source_complexes, indices)
    end
    return source_complexes
end
"""
sort_species(rn::ReactionSystem)

Sort species of a reaction network in an alphabetical order 
and return the indexes according to which they are sorted
"""
function sort_species(rn)
    if !(rn isa ReactionSystem)
        error("Input must be a Catalyst.ReactionSystem. Got: $(typeof(rn))")
    end
    sorted_species = sort(species(rn); by = x -> string(operation(x)))
    p = sortperm(species(rn); by = x -> string(operation(x)))
    return p
end

"""
ordering_complexes(rn::ReactionSystem)

Sort the complex of a reaction network in a given order, and return the incidence matrix sorted 
in consequence.
"""

function ordering_complexes(rn)
    S,Sub,Y,I,nb_complexes,nb_reactions,nb_species= extract_useful_from_crn(rn)
    rcs,B=reactioncomplexes(rn) #extract the reaction complexes
    p = sortperm(rcs) # sort the reactions complexes
    I = I[p,:] # return the reajusted incidence matrix, following this ordering.
    return I
end



"""
unique_vector_of_source_complexes(rn::ReactionSystem)

Returns a vector with each source complex of the given reaction network appearing once.
"""
function unique_vector_of_source_complexes(rn)
    symbolic_sources= [sum(s * x for (x, s) in zip(species(rn), col)) for col in eachcol(substoichmat(rn))]
    return unique(symbolic_sources)
end

"""
unique_vector_of_species(rn::ReactionSystem)

Returns vector of sorted species (in an alphabetical order)
"""
function unique_vector_of_species(rn)
    return sort(species(rn); by = x -> string(operation(x)))
end

"""
decide_identifiability(matrix::Matrix{Int64},source_complexes::Vector{Any})

Given a matrix (stoichiometric matrix or extended stoichiometric matrix) and source complexes,
 --> decide if the corresponding RN is identifiable
Actually it check linear independance of certain columns of the matrix.
Can be applied both to ODE and SDE
"""

function decide_identifiability(matrix, source_complexes)
    identifiable= true #start with identifiable =true
    extracted_columns= Matrix{Float64}[]
    for i in 1:length(source_complexes)
        # Extract columns de S corresponding to reactions where i is reactant
        push!(extracted_columns,matrix[:, source_complexes[i]])
        # if for one source complex the set of vector is linearly dependent then the ODE isn't identifiable
        if rank(extracted_columns[i]) != length(source_complexes[i])
            identifiable = false
        end 
    end   
    return identifiable
end

"""
solve_optimization_problem(matrix::Matrix{Int64},source_complexes::Vector{Any})

This function is applied when the reaction network isn't reaction-identifiable to find
two sets of rate constants giving the same ODEs/SDEs
It solves an optimization problem using Gurobi
"""

function solve_optimization_problem(matrice,source_complexes)
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, a[1:size(matrice,2)],Int, lower_bound=1) #first vector of rate constants, it should be integers
    @variable(model, b[1:size(matrice,2)],Int, lower_bound=1) #second vector of rate constants,it should be integers
    @constraint(model, sum((a[i] - b[i])^2 for i in 1:size(matrice,2))>=10e-1) # verify that the two vectors aren't the same
    for j in 1:length(source_complexes)
        if !isempty(source_complexes[j])
            #linear independence constraint
            @constraint(model, sum(a[i] * matrice[:,i] for i in source_complexes[j])== sum(b[i] * matrice[:,i] for i in source_complexes[j])) 
        end
    end
    #print(model)
    optimize!(model)
    return value.(a), value.(b)
end

"""
create_extended_matrix(rn::ReactionSystem,T::Matrix)

Create the matrix where each column is ((y'-y),(y'-y)(y'-y^T)).
It is used to to assess SDE identifiability, in fact we check linear independance of certain columns 
of this matrix, instead of linear independance of certain columns of the stoichiometric matrix.
"""

function create_extended_matrix(rn,T=nothing)
    S,Sub,Y,I,nb_complexes,nb_reactions,nb_species = extract_useful_from_crn(rn)
    if !(S isa AbstractMatrix{<:Real})
        error("Stoichiometry matrix 'S' must be a numeric matrix. Got: $(typeof(S))")
    end
    # create empty matrix of size nb_species + nb_species*nb_species,nb_reactions
    M = Matrix{Float64}(undef,nb_species + nb_species*nb_species,nb_reactions) 
    for j in 1:nb_reactions
        if isnothing(T)
            M[:,j] = vcat(S[:,j],reshape(S[:,j] * S[:,j]',(nb_species*nb_species,1)))
        else 
            M[:,j] = vcat(T*S[:,j],reshape(T*S[:,j] * S[:,j]'*T',(nb_species*nb_species,1)))
        end
    end
    return M
end

"""

symbolic_output(rn::ReactionSystem,a::Vector{Float64},b::Vector{Float64},equation_type::String,identifiable::Bool)

Function which returns symbolic output for identifiability where a and b are the rate constant vectors,
which gives the same ODEs or SDEs
"""

function symbolic_output(rn,a,b,equation_type::String,identifiable)
    S,Sub,Y,I,nb_complexes,nb_reactions,nb_species= extract_useful_from_crn(rn)
    Symbolics.@variables x[1:nb_species]
    x = [x[i] for i in 1:nb_species]
    vect =[]
    for i in 1:size(S,2)
        push!(vect,a[i]*prod(x.^(Sub[:,i])))
    end
    if equation_type == "ODE"
        if !identifiable
            println("This RN isn't reaction-identifiable w.r.t. its ODE : there exist 2 sets of rate constants k=",a, "and k'= ",b, " giving the same ODEs")
            display(latexstring("Let \$x_i\$ be the concentration of the \$i^\\text{th}\$ species. For both vectors of rate constants the drift vector is given by :  "))
            display(latexstring("A(x)= ",latexify(S*vect)))
        else 
            print("This RN is reaction-identifiable w.r.t. its ODE")
        end
    elseif equation_type == "SDE"
        if !identifiable
            println("The RN isn't reaction-identifiable w.r.t. its SDE. There exist two sets of rate constants giving same SDE : k=", value.(a), " and k'=", value.(b))
            display(latexstring("Let \$x_i\$ be the concentration of the \$i^\\text{th}\$ species. For both vectors of rate constants the drift vectors and diffusion matrices are given by :"))
            display(latexstring("A(x)= ",latexify(S*vect)))
            display(latexstring("B(x)= ",latexify(S*Diagonal(vect)*transpose(S))))
        else 
            print("This RN is reaction-identifiable w.r.t. its SDE")
        end
    else
        error("Invalid equation_type: '$equation_type'. Must be either 'ODE' or 'SDE'.")
    end
    return
end


""" 
create_rn_from_text

Function to create manually a reaction network 
"""

function create_rn_from_text(entree=nothing)
    if isnothing(entree)
        println("Creation of the RN : Enter reactions line by line")
        println("Just type reactions like 'A + B --> C + D'")
        println("Use '∅' for creation or degradation (e.g., '∅ --> A')")
        println("Enter an empty line to finish.")

        lignes = String[]
        seen_reactions = Set{String}()

        valid_species_regex = r"^(∅|\d*[A-Z]+(_\d+)?)$"

        norm = function(s)
            species = sort(split(strip(s), "+") .|> strip)
            for sp in species
                if !occursin(valid_species_regex, sp)
                    error("Invalid species name: '$sp'. Allowed: 'A', 'AX', 'A_1', or '∅'")
                end
            end
            return join(species, " + ")
        end

        counter = 1

        while true
            ligne = readline()
            ligne_clean = strip(ligne)
            isempty(ligne_clean) && break

            if !occursin(r"^\s*.*-->\s*.*$", ligne_clean)
                error("Invalid format: '$ligne_clean'. Expected: 'A + B --> C + D'")
            end

            rp_split = split(ligne_clean, "-->")
            if length(rp_split) != 2
                error("Reaction must contain exactly one '-->'")
            end
            reactants, products = rp_split .|> strip

            normalized_reaction = norm(reactants) * " --> " * norm(products)

            if normalized_reaction in seen_reactions
                error("Duplicate reaction detected: '$normalized_reaction'")
            end

            rate = "k_$(counter)"
            push!(seen_reactions, normalized_reaction)
            push!(lignes, "$rate, $normalized_reaction")
            counter += 1
        end

        entree = join(lignes, "\n    ")
    end

    model_code = """
    @reaction_network begin
        $entree
    end
    """
    
    return eval(Meta.parse(model_code))
end

"""
ode_reaction_identifiability(rn::ReactionSystem)

Output whether a reaction network is identifiable w.r.t. its ODE or not.
If not provide two vectors of rate constants giving the same ODE and the corresponding drift vector.
"""

function ode_reaction_identifiability(rn=nothing)

    #if no input, let the user enter manually his reaction network
    if isnothing(rn)
        rn = create_rn_from_text()
    end
    #verify type of input
    if !(rn isa Catalyst.ReactionSystem)
        error("Expected a Catalyst.ReactionSystem. Got: $(typeof(rn))")
    end
    #display reaction network
    println( "The reaction network is : ")
    display(latexify(rn))
    numreactions(rn) > 0 ||
        error("There must be at least one reaction.")
    S,Sub,Y,I,nb_complexes,nb_reactions,nb_species= extract_useful_from_crn(rn)
    #define a vector where will be stocked source complexes
    source_complexes = find_source_complexes(rn)
    identifiable = decide_identifiability(S,source_complexes)
    #if identifiable == false, solve the optimization problem and find 2 vectors of rate constants giving same ODEs
    if identifiable == false
        a,b = solve_optimization_problem(S,source_complexes)
        symbolic_output(rn,a,b,"ODE",identifiable)
    else
        print("This RN is reaction-identifiable w.r.t. its ODE")
    end
    return 
end

"""
sde_reaction_identifiability(rn::ReactionSystem)

Output whether a reaction network is identifiable w.r.t. its SDE or not.
If not provide two vectors of rate constants giving the same generators and the corresponding drift vector and diffusion matrix.
"""

function sde_reaction_identifiability(rn=nothing)

    #if no input, let the user enter manually his reaction network
    if isnothing(rn)
        rn = create_rn_from_text()
    end
    #verify type of input
    if !(rn isa Catalyst.ReactionSystem)
        error("Expected a Catalyst.ReactionSystem. Got: $(typeof(rn))")
    end
    #display reaction network
    println( "The reaction network is : ")
    display(latexify(rn))
    numreactions(rn) > 0 ||
        error("There must be at least one reaction.")
    S,Sub,Y,I,nb_complexes,nb_reactions,nb_species= extract_useful_from_crn(rn)
    #define a vector where will be stocked source complexes
    source_complexes = find_source_complexes(rn)
    M = create_extended_matrix(rn)
    identifiable = decide_identifiability(M,source_complexes)
    #if identifiable == false, solve the optimization problem and find 2 vectors of rate constants giving same generators.
    if identifiable == false
        a,b = solve_optimization_problem(M,source_complexes)
        symbolic_output(rn,a,b,"SDE",identifiable)
    else
        print("This RN is reaction-identifiable w.r.t. its SDE")
    end
    return
end


"""
solve_confoundability_optimization_problem(rn1::ReactionSystem,rn_2::ReactionSystem,equation_type::String,conjugacy::Bool)

Solve optimization problem for confoundablity or linear conjugacy both w.r.t. their ODEs or SDEs
Print directly the results.
"""


function solve_confoundability_optimization_problem(rn_1,rn_2,equation_type::String,conjugacy::Bool)
    #verify that the source complexes are the same, if not it isn't confoundable
    source_complexes_1 = find_source_complexes(rn_1)
    source_complexes_2 = find_source_complexes(rn_2)
    source_complexes_1 = filter(x -> !isempty(x), source_complexes_1)
    source_complexes_2 = filter(x -> !isempty(x), source_complexes_2)
    #have_same_source_complexes(source_complexes_1,source_complexes_2)
    S_1,Sub_1,Y_1,I_1,nb_complexes_1,nb_reactions_1,nb_species_1= extract_useful_from_crn(rn_1)
    S_2,Sub_2,Y_2,I_2,nb_complexes_2,nb_reactions_2,nb_species_2= extract_useful_from_crn(rn_2)
    model = Model(Gurobi.Optimizer)
    S_1 = S_1[sort_species(rn_1),:]
    S_2 = S_2[sort_species(rn_2),:]
    #see what to do with complexes ordering
    # I_1 = ordering_complexes(rn_1)
    # I_2 = ordering_complexes(rn_2)
    set_silent(model)
    set_optimizer_attribute(model, "TimeLimit", 60)  # limit to 60 seconds
    set_optimizer_attribute(model, "MIPGap", 0.05)   # allow 5% suboptimality
    set_optimizer_attribute(model, "NonConvex", 2)   # needed for c[1]*c[2]
    #if we solve a linear conjugacy problem, we add some variables c of the diagonal matrix
    if conjugacy 
        @variable(model,c[1:size(S_2,1)])
        for i in 1:size(S_2,1)
            @constraint(model,c[i] >= 10e-3)
        end
    end
    # if SDE, we don't work with columns of the stoichiometric matrix, but with ones of the extended stoichiomatric matrix
    if equation_type == "SDE"
        # for first RN, it is simply the extended stoichiomatric matrix
        M_1 = create_extended_matrix(rn_1)
        if conjugacy 
            # if we look for linear conjugacy then for second RN, the diagonal matrix c comes into action
            M_2 = Matrix{Any}(undef,nb_species_2 + nb_species_2*nb_species_2,nb_reactions_2) 
            for j in 1:nb_reactions_2
                M_2[:,j] = vcat(diagm(c)*S_2[:,j],reshape(diagm(c)*S_2[:,j] * S_2[:,j]'*diagm(c)',(nb_species_2*nb_species_2,1)))
            end  
            #M_2 = map(safe_simplify_expr,M_2)
        else 
            # if we look for confoundablity then for second RN, it is simply the extended stoichiomatric matrix
            M_2 = create_extended_matrix(rn_2)
        end
        matrice_1 = M_1
        matrice_2 = M_2
    elseif equation_type=="ODE"
        #same things for ODE
        matrice_1 = S_1
        if conjugacy 
            matrice_2 = diagm(c)*S_2
        else
            matrice_2 = S_2
        end
    else
        error("Invalid equation_type: '$equation_type'. Must be either 'ODE' or 'SDE'.")
    end
    # rate constants constants must be integers and not zeros
    @variable(model, a[1:size(matrice_1,2)],Int, lower_bound=1)
    @variable(model, b[1:size(matrice_2,2)],Int, lower_bound=1)
    @constraint(model, sum((a[i])^2 for i in 1:size(matrice_1,2))>=10e-1)
    @constraint(model, sum((b[i])^2 for i in 1:size(matrice_2,2))>=10e-1)
    for j in 1:length(source_complexes_1)
        if !isempty(source_complexes_1[j])
            left_exprs = sum(a[i] * matrice_1[:, i] for i in source_complexes_1[j])
            right_exprs = sum(b[i] * matrice_2[:, i] for i in source_complexes_2[j])
            @constraint(model, left_exprs .== right_exprs)
        end
    end
    #print(model)
    optimize!(model)
    if termination_status(model)==MOI.OPTIMAL
        Symbolics.@variables x[1:nb_species_1]
        x = [x[i] for i in 1:nb_species_1]
        vect_1 =[]
        for i in 1:size(S_1,2)
            push!(vect_1,value.(a)[i]*prod(x.^(Sub_1[:,i])))
        end
        vect_2 =[]
        for i in 1:size(S_2,2)
            push!(vect_2,value.(b)[i]*prod(x.^(Sub_2[:,i])))
        end
        if equation_type == "ODE"         
            if conjugacy 
                println("The 2 RNs are linearly-conjugated w.r.t. their ODEs")
                println("The two RNs have same ODEs for rate constants : k=",value.(a), "and k'=",value.(b), " and c = ", value.(c))
            else
                println("The 2 RNs are confoundable w.r.t. their ODEs")
                println("The two RNs have same ODEs for rate constants : k=",value.(a), "and k'=",value.(b))
            end
            display(latexstring("Let \$x_i\$ be the concentration of your \$i^\\text{th}\$ first species. For both vectors of rate constants the drift vectors is given by :  "))
            display(latexstring("A(x) = ",latexify(S_1*vect_1)))
        else        
            if conjugacy 
                println("The 2 RNs are linearly-conjugated w.r.t. their SDEs")
                println("The two RNs have same generator for rates : k=",value.(a), "and k'=",value.(b), " and c = ", value.(c))
            else 
                println("The 2 RNs are confoundable w.r.t. their SDEs")
                println("The two RNs have same generator for rates : k=",value.(a), "and k'=",value.(b))
            end
                       
            display(latexstring("Let \$x_i\$ be the concentration of your \$i^\\text{th}\$ first species. For both vectors of rate constants the drift vectors and diffusion matrices are given by :  "))
            display(latexstring("A(x) = ",latexify(S_1*vect_1)))
            display(latexstring("B(x) = ",latexify(S_1*Diagonal(vect_1)*transpose(S_1))))
        end
    elseif termination_status == MOI.TIME_LIMIT
        println("Optimization stopped due to time limit.There is no intersection between the two cones. The 2 RNs are not confoundable w.r.t. their SDEs")
    else 
        if equation_type == "ODE"
            println("The 2 RNs are not confoundable w.r.t. their ODEs")
        else
            println("The 2 RNs are not confoundable w.r.t. their SDEs")
        end
    end
end

"""
ode_confoundability(rn1::ReactionSystem,rn_2::ReactionSystem)

Output if two reactions networks are confoundable w.r.t. their ODEs or not.
If so it gives two sets of rate constants vectors giving same ODEs
"""


function ode_confoundability(rn_1=nothing, rn_2=nothing)

    #if no input, let the user enter manually a RN
    if isnothing(rn_1)
        rn_1 = create_rn_from_text()
    end
    #display first RN
    println( "First reaction network : ")
    display(latexify(rn_1))
    #if no input, let the user enter manually a RN
    if isnothing(rn_2)
        rn_2 = create_rn_from_text()
    end
    #display second RN
    println( "Second reaction network : ")
    display(latexify(rn_2))
    numreactions(rn_1) > 0 ||
        error("There must be at least one reaction.")
    numreactions(rn_2) > 0 ||
        error("There must be at least one reaction")
    if Set(unique_vector_of_source_complexes(rn_1)) !=Set(unique_vector_of_source_complexes(rn_2)) 
        error("The two RNs don't have same source complexes")
    end
    if Set(unique_vector_of_species(rn_1)) != Set(unique_vector_of_species(rn_2))
        error("The two RNs don't have same source species")
    end
    solve_confoundability_optimization_problem(rn_1,rn_2,"ODE",false)
    return
end


"""
ode_linear_conjugacy(rn1::ReactionSystem,rn_2::ReactionSystem)

Output if two reactions networks are linearly conjugated w.r.t. their ODEs or not.
If so it gives two sets of rate constants vectors, and of scaling constants, giving same ODEs.
"""

function ode_linear_conjugacy(rn_1=nothing, rn_2=nothing)

    if isnothing(rn_1)
        rn_1 = create_rn_from_text()
    end
    println( "First reaction network : ")
    display(latexify(rn_1))
    if isnothing(rn_2)
        rn_2 = create_rn_from_text()
    end
    println( "Second reaction network : ")
    display(latexify(rn_2))
    numreactions(rn_1) > 0 ||
        error("There must be at least one reaction.")
    numreactions(rn_2) > 0 ||
        error("There must be at least one reaction.")
    if Set(unique_vector_of_source_complexes(rn_1)) !=Set(unique_vector_of_source_complexes(rn_2)) 
        error("The two RNs don't have same source complexes")
    end
    if Set(unique_vector_of_species(rn_1)) != Set(unique_vector_of_species(rn_2))
        error("The two RNs don't have same species")
    end
    solve_confoundability_optimization_problem(rn_1,rn_2,"ODE",true)
    return
end

"""
sde_confoundability(rn1::ReactionSystem,rn_2::ReactionSystem)

Output if two reactions networks are confoundable w.r.t. their SDEs or not.
If so it gives two sets of rate constants vectors giving same generators.
"""

function sde_confoundability(rn_1=nothing, rn_2=nothing)
    if isnothing(rn_1)
        rn_1 = create_rn_from_text()
    end
    println( "First reaction network : ")
    display(latexify(rn_1))
    if isnothing(rn_2)
        rn_2 = create_rn_from_text()
    end
    println( "Second reaction network : ")
    display(latexify(rn_2))
    numreactions(rn_1) > 0 ||
        error("There must be at least one reaction.")
    numreactions(rn_2) > 0 ||
        error("There must be at least one reaction.")   

    
    if Set(unique_vector_of_source_complexes(rn_1)) !=Set(unique_vector_of_source_complexes(rn_2)) 
        error("The two RNs don't have same source complexes")
    end
    if Set(unique_vector_of_species(rn_1)) != Set(unique_vector_of_species(rn_2))
        error("The two RNs don't have same species")
    end
    solve_confoundability_optimization_problem(rn_1,rn_2,"SDE",false)
    return
end

"""
sde_linear_conjugacy(rn1::ReactionSystem,rn_2::ReactionSystem)

Output if two reactions networks are linearly conjugated w.r.t. their SDEs or not.
If so it gives two sets of rate constants vectors, and of scaling constants, giving same generators.
"""

function sde_linear_conjugacy(rn_1=nothing, rn_2=nothing)

    if isnothing(rn_1)
        rn_1 = create_rn_from_text()
    end
    println( "First reaction network : ")
    display(latexify(rn_1))
    if isnothing(rn_2)
        rn_2 = create_rn_from_text()
    end
    println( "Second reaction network : ")
    display(latexify(rn_2))
    numreactions(rn_1) > 0 ||
        error("There must be at least one reaction.")
    numreactions(rn_2) > 0 ||
        error("There must be at least one reaction.")
    if Set(unique_vector_of_source_complexes(rn_1)) !=Set(unique_vector_of_source_complexes(rn_2)) 
        error("The two RNs don't have same source complexes")
    end
    if Set(unique_vector_of_species(rn_1)) != Set(unique_vector_of_species(rn_2))
        error("The two RNs don't have same species")
    end
    solve_confoundability_optimization_problem(rn_1,rn_2,"SDE",true)
    return
end

