using DataStructures # for queue


mutable struct SupportedSolutionPoint
    sp :: VectorSolutionPoint
    λ :: Vector{Float64}
    is_integer :: Bool 

end



# ----------------------------------
# ---------- Node ------------------
# ----------------------------------

mutable struct Node 
    num::Int64                    
    depth::Int64                # depth in tree
    pred::Node                  # predecessor
    succs::Vector{Node}         # successors
    var_idx::Int64                  # index of the chosen variable to be split
    var_bound::Int64            # variable bound
    bound_type :: Int64         # 1 : <= ; 2 : >= 
    activated::Bool             # if the node is active
    pruned::Bool                # if the node is pruned
    pruned_type::PrunedType      # if the node is fathomed, restore pruned type
    deleted::Bool               # if the node is supposed to be deleted
    lower_bound_set::Vector{SupportedSolutionPoint}        # local lower bound set    
    assignment::Vector{Tuple{Int64, Int64, Int64}}  # (varidex, varbound, boundtype)

    Node() = new()

    function Node(num::Int64, depth::Int64 ;
        pred::Node=Node(), succs::Vector{Node}=Vector{Node}(), var_idx::Int64=0, var_bound::Int64=0, bound_type::Int64=0
   )
        n = new()
        n.num = num
        n.depth = depth
        n.pred = pred
        n.succs = succs
        n.var_idx = var_idx
        n.var_bound = var_bound
        n.bound_type = bound_type

        n.activated = true 
        n.pruned = false
        n.pruned_type = NONE
        n.deleted = false
        n.lower_bound_set = Vector{SupportedSolutionPoint}()
        n.assignment = Vector{Tuple{Int64, Int64, Int64}}()

        f(t) = nothing 
        finalizer(f, n)
    end
end

function Base.:show(io::IO, n::Node)
    println(io, "\n\n # ----------- node $(n.num) : \n", 
    "depth = $(n.depth) \n",
    "pred = $(n.pred.num) \n",
    "var[ $(n.var_idx) ] = $(n.var_bound) \n",
    "activated = $(n.activated) \n",
    "pruned = $(n.pruned) \n",
    "prunedType = $(n.prunedType)"
    )
    print(io, "succs = [ ")
    for s in n.succs print(io, "$(s.num), ") end
    println(io, " ]")

    print(io, "LBS = [ ")
    for s in n.lower_bound_set
        print(io, "$(s.sp.y) , ")
    end
    println(io, "] ")
end


"""
Return `true` if the given node is the root of a branch-and-bound tree.
"""
function isRoot(node::Node)
    return node.depth == 0 # !isdefined(node, :pred) 
end

"""
Return `true` if the given `node` has activated/non-explored successor(s).
"""
function hasNonExploredChild(node::Node)
    for c in node.succs
        if c.activated return true end 
    end
    return false
end

"""
Delete the given node in B&B tree. (should be a private function)
"""
function Base.delete!(node::Node)
    node.deleted = true ; node = nothing               # remove from the memory
end

"""
Prune the given node in a B&B tree and delete all successors of the pruned node.
"""
function prune!(node::Node, reason::PrunedType)
    node.pruned = true
    node.pruned_type = reason

    to_delete = node.succs[:]
    node.succs = Vector{Node}()

    while length(to_delete) > 0
        n = pop!(to_delete)
        to_delete = vcat(to_delete, n.succs[:])
        delete!(n)
    end
end


"""
From the actual node, go up to the root to get the partial assignment of variables.
"""
function getPartialAssign(actual::Node)::Vector{Tuple{Int64, Int64, Int64}}
    assignment = Vector{Tuple{Int64, Int64, Int64}}() # (varidex, varbound, boundtype)
    if isRoot(actual) # the actual node is the root 
        return assignment
    end
    predecessor = actual.pred
    push!(assignment, (actual.var_idx, actual.var_bound, actual.bound_type))
    # todo 
    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        if !actual.EPB assignment[actual.var] = actual.var_bound end 
    end
    return assignment
end


# ----------------------------------
# ---------- branching -------------
# ----------------------------------
