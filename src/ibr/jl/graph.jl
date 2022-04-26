######################################################################
# graph algorithms
######################################################################

"""
    A graph is represented by:

    - An edge list `fel`: a table with columns "source vertex",
      "destination vertex"

    - A backward edge list `bel`: a table with a single column
      "forward edge number"

    - A "forward star" list `fs` with the start indices of the forward
      egdes for each vertex

    - A "backward star" list `bs` with the start indices of the backward
      egdes for each vertex
"""
struct Graph
    ### Note: DS is Ok for directed graphs, most alg below assume undirected
    
    # number of vertices and edges
    n :: Int ## @assert n == max(maximum(fel[:,1]),maximum(fel[:,2]))
    m :: Int ## @assert m == size(fel,1)
    # forward edge list (src,dst)
    fel :: Array{Int,2}
    # forward star
    fs :: Array{Int}
    # backward edge list (index into fel)
    bel :: Array{Int}
    # backward star
    bs :: Array{Int}
    ## some helpers: scratch memory for algorithm
    mark :: BitArray{1}
    idx :: Vector{Int}
end

## DIMACS format
function readInstance(io::IO)
    ## (1) read the file
    n = m = 0
    e = 0
              
    while !eof(io)
        l = readline(io)
        occursin(r"^c ",l) && continue
        ma = match(r"^p\s+edge\s+(\d+)\s+(\d+)",l)
        if ma != nothing
              n = parse(Int,ma[1])
              m = parse(Int,ma[2])
              fel = zeros(Int,2*m,2)
        end
        ma = match(r"^e\s+(\d+)\s+(\d+)",l)
        if ma != nothing
            e += 1
            fel[e,1] = parse(Int,ma[1])
            fel[e,2] = parse(Int,ma[2])
            e += 1
            fel[e,1] = fel[e-1,2]
            fel[e,2] = fel[e-1,1]
        end
    end
    fel=fel[1:e,:]
    @assert e == 2*m || e == m
    
    # (2) create forward star pointers
    fel = unique(sortslices(fel,dims=1),dims=1)
    @assert size(fel,1) % 2 == 0
    m = size(fel,1)/2
    fs=map(x->findfirst(isequal(x),fel[:,1]),1:n+1)
    fs[n+1]=2*m+1
    for i=n:-1:1
        if fs[i]==0
            fs[i]=fs[i+1]
        end
    end

    # (3) create backward star pointers
    bel=hcat(fel,1:2*m) ## add forward edge identifier
    bel=sortslices(bel,dims=1,lt=(a,b)->a[2]<b[2]) ## sort by (destination)
    bs=map(x->findfirst(isequal(x),bel[:,2]),1:n+1)
    bs[n+1]=2*m+1
    for i=n:-1:1
        if bs[i]==0
            bs[i]=bs[i+1]
        end
    end
    # (4.1) drop all columns, except pointer to corresponding forward edge
    bel=bel[:,3]

    Graph(n,m,fel,fs,bel,bs,falses(n),zeros(Int,n))
end

## apply function to edges in subgraph induced by vertices V
function applyEdges(G::Graph,V::Vector{Int},f)
    fill!(G.mark,false)
    G.mark[V] .= true
    #println("Vertices $V")
    for v in V
        for oa=G.fs[v]:G.fs[v+1]-1
            u=G.fel[oa,2]
            (G.mark[u] && u<v) || continue ## check later: v<u better?
            f(v,u)
        end
    end
end

## count edges in subgraph induced by vertices V
function edges(G::Graph,V::Vector{Int})
    ie = 0
    fill!(G.mark,false)
    G.mark[V] .= true
    for v in V
        for oa=G.fs[v]:G.fs[v+1]-1
            u=G.fel[oa,2]
            (G.mark[u] && v<u) || continue
            ie += 1
        end
    end
    ie
end

### find the number of reachable vertices from `v0` (by BFS)
function bfs(G::Graph,v0;depth=zeros(Int,G.n))
    visited=zeros(Int,G.n) ## 0: not, 1: queue, 2: done
    
    Q = [v0]
    bipartite = true

    r = 0
    while length(Q)>0
        v = popfirst!(Q)
        visited[v]=2
        r += 1
        for oa=G.fs[v]:G.fs[v+1]-1
            w = G.fel[oa,2]
            if visited[w]==0
                visited[w]=1
                depth[w]=depth[v]+1
                push!(Q,w)
            else
                if depth[w]==depth[v]
                    bipartite = false
                end
            end
        end
    end
    (r,bipartite,findall(x->x!=0,visited))
end

function isBipartite(G::Graph)
    bfs(G::Graph,1)[2]
end

function maxDeg(G::Graph)
    md = 0
    for v in 1:G.n
        md = max(md,G.fs[v+1]-G.fs[v])
    end
    md
end

function isOddCycle(G::Graph)
    G.n != G.m && return false
    for v in 1:G.n
        G.fs[v+1]-G.fs[v] != 2 && return false
    end
    bfs(G,1)[1] == G.n
end

function isComplete(G::Graph)
    2*G.m == G.n*(G.n-1)
end

## bounds on chromatic number (somewhat inefficient, calls bfs possibly twice)
function boundχ(G::Graph)
    G.n==0 && return (0,0)
    G.n==1 && return (1,1)
    G.m==0 && return (1,1)
    isBipartite(G) && return (2,2)
    isOddCycle(G) && return (3,3)
    isComplete(G) && return (G.n,G.n)
    (3,maxDeg(G))  ## trivial, Brook's theorem
end

## write graph induced by a partition in DOT format
function writeDOT(io::IO,G::Graph,P::Partition)
    println(io,"graph {")
    V=find(P.crit.|P.open)
    for v in V
        color=P.crit[v] ? "red" : "yellow"
        println(io,"$v [color=$(color)]")
    end
    applyEdges(G,V,(u,v)->println(io,"$v -- $u"))
    println(io,"}")
end

## write graph induced by vertex set V in DIMACS format
function writeDIMACS(io::IO,G::Graph,V::Vector{Int};comment::String="")
    fill!(G.idx,0)
    cv = 0
    comment != "" && println(io,"c Comment: $comment")
    comment != "" && print(io,"c vertices ")
    for v in sort(V)
        comment != "" && print(io,"$v ")
        cv += 1
        G.idx[v]=cv
    end
    comment != "" && println(io)
    ne = edges(G,V)
    println(io,"p edge $(length(V)) $ne")
    applyEdges(G,V,(u,v)->println(io,"e $(G.idx[v]) $(G.idx[u])"))
end
