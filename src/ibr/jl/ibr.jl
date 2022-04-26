#!/usr/bin/env julia
######################################################################
# a draft implementation in Julia
######################################################################
using ArgParse
using DataStructures
using Printf
using DelimitedFiles

# partition into crit, open, disp
struct Partition
    ## vertex status
    crit :: BitArray{1}
    open :: BitArray{1}
    # @assert crit .&& open == false
    # disp = !crit && !open

    Open :: PriorityQueue{Int,Tuple{Int,Int}}
    Disp :: Stack{Int}

    f :: Vector{Int} # flips

    function Partition(G)
        Q = PriorityQueue{Int,Tuple{Int,Int}}()
        for v=1:G.n
            enqueue!(Q,v,(0,G.fs[v+1]-G.fs[v]))
        end
        new(falses(G.n),trues(G.n),Q,Stack{Int}(),zeros(Int,G.n))
    end
end

include("graph.jl")

@enum Algorithms TabuCol HybridEA PCTC

## color using TabuCol, return number of colors `x` and the coloring
function color_tabuCol(G::Graph,V::Vector{Int})
    cfile="tabuCol-$(getpid()).sol"
    pr=open(`$(base)/graphCol/color --maxIt 20000 --coloring $(cfile) /proc/self/fd/0`,"r+")
    writeDIMACS(pr,G,V)
    close(pr.in)
    x=parse(Int,readline(pr))
    close(pr)
    r=(x,readdlm(cfile,Int,skipstart=1)[:,2].+1)
    rm(cfile)
    r
end

## color using HybridEA
function color_hybridea(G::Graph,V::Vector{Int},k::Int)
    cfile="hybridea-$(getpid()).sol"
    pr=open(`$(base)/gCol/HybridEA/HybridEA --target $k --printcolors --solution $(cfile) /proc/self/fd/0`,"r+")
    writeDIMACS(pr,G,V)
    close(pr.in)
    x=parse(Int,readline(pr))
    close(pr)
    r=(x,readdlm(cfile,Int,skipstart=1)[:,1].+1)
    rm(cfile)
    r
end

## color using HybridEA
function color_pctc(G::Graph,V::Vector{Int},k::Int)
    cfile="pctc-$(getpid()).sol"
    pr=open(`$(base)/gCol/PartialColAndTabuCol/PartialColAndTabuCol --target $k --printcolors --solution $(cfile) /proc/self/fd/0`,"r+")
    writeDIMACS(pr,G,V)
    close(pr.in)
    x=parse(Int,readline(so))
    close(pr)
    r=(x,readdlm(cfile,Int,skipstart=1)[:,1].+1)
    rm(cfile)
    r
end

## color graph induced by vertex set V
function color(G::Graph,V::Vector{Int},k::Int;algorithm=TabuCol)
    ## (1) remove some bases cases
    length(V)==0 && return 0
    length(V)==1 && return 1
    ne = edges(G,V)
    ne == 0 && return 1
    length(V)==2 && return ne+1
    ## (2) check easy targets
    k==0 && return G.n
    k==1 && return G.n
    k==2 && isBipartite(G) ? 2 : G.n

    ## (3) external coloring
    if algorithm==TabuCol
        color_tabuCol(G,V)[1]
    elseif algorithm==HybridEA
        color_hybridea(G,V,k)[1]
    else
        @assert algorithm==PCTC
        color_pctc(G,V,k)[1]
    end
end

## recolor after removal of v from V
function recolor_down(G::Graph,V::Vector{Int},v::Int,k::Int;algorithm=TabuCol)
    color(G,V,k,algorithm=algorithm)
end

## recolor after inclusion of v from V
function recolor_up(G::Graph,V::Vector{Int},v::Int,k::Int;algorithm=TabuCol)
    color(G,V,k,algorithm=algorithm)
end

const disp_open = (0,1)
const disp_crit = (1,0)
const open_disp = .-disp_open
const open_crit = (1,-1)
const crit_disp = .-disp_crit
const crit_open = .-open_crit

function moveCritOpen(G::Graph,P::Partition,v::Int)
    P.crit[v]=false
    P.open[v]=true
    prio=updatePrioritiesGet(G,P,v,crit_open)
    enqueue!(P.Open,v,prio)
    P.f[v]+=1
end

function nextOpenDisp(G::Graph,P::Partition)
    v=dequeue!(P.Open)
    P.open[v]=false
    push!(P.Disp,v)
    updatePriorities(G,P,v,open_disp)
    P.f[v]+=1
    v
end

function nextDispOpen(G::Graph,P::Partition)
    v=pop!(P.Disp)
    P.open[v]=true
    prio=updatePrioritiesGet(G,P,v,disp_open)
    enqueue!(P.Open,v,prio)
    P.f[v]+=1
    v
end

function nextDispCrit(G::Graph,P::Partition)
    v = pop!(P.Disp)
    P.crit[v]=true
    updatePriorities(G,P,v,disp_crit)
    P.f[v] += 1
    v
end

## v: vertex that changes the group
## delta: change vector, see constants above
## return: priority of v
function updatePrioritiesGet(G::Graph,P::Partition,v::Int,delta::Tuple{Int,Int})
    prio=(0,0)
    for oa=G.fs[v]:G.fs[v+1]-1
        u=G.fel[oa,2]
        if P.crit[u] prio = prio .+ (1,0) end
        if P.open[u]
            prio = prio .+ (0,1)
            P.Open[u] = P.Open[u] .+ delta
        end
    end
    prio
end

## as above, but don't need priority of v
function updatePriorities(G::Graph,P::Partition,v::Int,delta::Tuple{Int,Int})
    for oa=G.fs[v]:G.fs[v+1]-1
        u=G.fel[oa,2]
        if P.open[u]
            P.Open[u] = P.Open[u] .+ delta
        end
    end
end

function ibr(G::Graph,k::Int,R::Int,p::Float64,ηf::Float64;verbose::Bool=false)
    ## (1) prepare data structures
    P = Partition(G)
    Crit = Vector{Int}[]
    e = zeros(Int,G.n) # elite occurences

    ## (2) repeatedly search for a k-critical component
    cc = 0
    for r = 1:R
        verbose && println("===== Round $r =====")

        ## (2.1) find a k-critical set
        verbose && println("===== Calling findkVCS $r =====")
        success = findkVCS(G,P,cc,k,verbose=verbose)
        !success && return Int[]

        ## (2.2) add to "elite" set
        push!(Crit,findall(P.crit))
        e[Crit[end]] .+= 1

        ## smallest graph with χ=k is complete, so if complete, we may stop
        ncrit = length(findall(P.crit))
        mcrit = edges(G,findall(P.crit))
        verbose && println("Elite: $ncrit $mcrit")
        ncrit == k && ncrit*(ncrit-1) == 2*mcrit && break

        ## (2.3) perturb and update
        η = ceil(Int,length(Crit[end])*ηf)
        S = perturb(G,P,Crit,e,η,p)
        for v in S
            moveCritOpen(G,P,v)
        end
        cc = color(G,findall(P.crit),k-1)
        if cc<k
            # move all disp to open
            while !isempty(P.Disp)
                nextDispOpen(G,P)
            end
        else
            ## TBD: have better critical graph here: record it!
            # move all open to disp
            while !isempty(P.Open)
                nextOpenDisp(G,P)
            end
            # move all crit to open
            for v in findall(P.crit)
                moveCritOpen(G,P,v)
            end
            cc = 0
        end
    end

    ## find smallest critical subgraph
    bi = 0
    b = (typemax(Int),typemax(Int))
    for r = 1:R
        c=(length(Crit[r]),edges(G,Crit[r]))
        if c<b
            b=c
            bi=r
        end
    end
    Crit[bi]
end

function findkVCS(G::Graph,P::Partition,cc::Int,k::Int; verbose::Bool=false)
    verbose && println("===== findkVCS")
    verbose && println("     Start    : crit=$(length(findall(P.crit)))/$cc open=$(length(P.Open))/?")
    img=1
    while cc<k
        co = cc
        if !isempty(P.Open)
            #println("Top prio: $(peek(P.Open))")
            v=nextOpenDisp(G,P)
            co = recolor_down(G,findall(P.crit .| P.open),v,k-1,algorithm=TabuCol)
	    verbose && println("     Remove $v: crit=$(length(findall(P.crit)))/$cc open=$(length(P.Open))/$co")
        end
        ## Else: co has been judged to be >=k-critical, but cc not (in line (RU)), even though they are the
        ## same. The default "co = cc" from above applies.
        first_iteration = true
        w = 0
        while co < k
            if first_iteration
                first_iteration = false
                w = nextDispCrit(G,P)
            else
                isempty(P.Disp) && return false  ## no example
                nextDispOpen(G,P)
            end
            # here we color the original open+crit again; but since we're run a heuristic the # of colors
            # could drop, so we call the heuristic again. A second iteration of the while loop is only needed
            # if it indeed does drop.  TBD: if, however, the number of colors increases, we should use the
            # former.
            co = recolor_up(G,findall(P.crit .| P.open),w,k-1,algorithm=TabuCol)
            verbose && println("     Add    $w: crit=$(length(findall(P.crit)))/$cc open=$(length(P.Open))/$co")
        end

        # if crit did not increase, we have a witness that it is <k-colorable, i.e. cc<k, so no need to recolor
        if w>0
            cc = recolor_up(G,findall(P.crit),w,k-1,algorithm=HybridEA) ## (RU)
            verbose && println("     Crit   $w: crit=$(length(findall(P.crit)))/$cc open=$(length(P.Open))/$co")
        end
        # Note: if v is the only vertex of its color, then it is obviously critical, and we can move it directly and proceed to newcrit
        #if length(findall(P.crit))+length(P.Open)<25 && verbose
        #    println("Writing image $(img).")
        #    open("img$(img).dot","w") do io
        #        writeDOT(io,G,P)
        #    end
        #    img+=1
        #end
    end
    @assert length(findall(P.crit))>=k
    true
end

# select critical vertices to move back to open from last critical set
function perturb(G::Graph,P::Partition,Crit,e::Vector{Int},η::Int,p::Float64)
    ## (1) compute scores
    fmax = maximum(P.f)
    score = zeros(G.n)
    for v in Crit[end]
        score[v] = e[v]/length(Crit)*(1-e[v]/length(Crit)) + 1-P.f[v]/fmax
    end
    ## (2) select vertices
    π = sortperm(Crit[end],lt=(a,b)->score[a]>score[b])
    @assert η>0
    jsum = sum(1/j for j=1:η)
    cj = 0
    for j in 1:η
        if 1>j*rand() # j*jsum*rand()
            cj += 1
            π[cj]=π[j]
        end
    end
    Crit[end][π[1:cj]]
end


# (0) check the input
s = ArgParseSettings()
@add_arg_table s begin
    "instance"
    help = "instance"
    required = true
    "--k"
    help = "number of critical vertices k, default 0 means search by increasing values of k"
    arg_type = Int
    default = 0
    "--R"
    help = "number of repetitions"
    arg_type = Int
    default = 20
    "--p"
    help = "acceptance probability after perturbation"
    arg_type = Float64
    default = 0.1
    "--etaf"
    help = "fraction of critical vertices considered during perturbation"
    arg_type = Float64
    default = 0.5
    "--solution"
    help = "solution output file, default ibr.sol"
    default = "ibr.sol"
    "--base"
    help = "base directory for external programs"
    default = "/home/mrpritt/Home/ibr"
    "--verbose"
    help = "show verbose output"
    action = :store_true
end

args = parse_args(s)

k=args["k"]
name=basename(args["instance"])
base=args["base"]

## (1) read the instance, prepare
G = open(args["instance"],"r") do io
    readInstance(io)
end
(lχ,uχ)=boundχ(G)
println("# bounds on χ: $lχ $uχ")

## (2) find k-critical subgraphs (vertex list V)
start=time()
V = Int[]
if k>0
    V = ibr(G,k,args["R"],args["p"],args["etaf"],verbose=args["verbose"])
    @printf("%s %3d %4d %5d %4d %5d %5.1f\n",name,k,G.n,G.m,length(V),edges(G,V),time()-start)
else
    k=1
    for k=1:G.n
        iterstart = time()
        V′ = ibr(G,k,args["R"],args["p"],args["etaf"],verbose=args["verbose"])
        length(V′)==0 && break
        V=V′
        @printf("%s %3d %4d %5d %4d %5d %5.1f\n",name,k,G.n,G.m,length(V),edges(G,V),time()-iterstart)
    end
end

## (3) print solution
open(args["solution"],"w") do io
    writeDIMACS(io,G,V,comment="source graph $name")
end
