using Oscar

"""
    has_Mukai_edge_graph(M::Matrix{Int64},deg_mon::Vector{Vector{Int64}})

    The inputs are the 129x14 matrix M of the exponents of the initial forms of the generators of the Cox ring and a vector deg_mon of multidegress of the generators.
    The output is a boolean which is true if and only if the polytope obtained as the convex hull of the rows of M has the Mukai edge-graph.

"""

function has_Mukai_edge_graph(M,deg_mon)
    P=convex_hull(M);
    EdGr=Oscar.Graphs.edgegraph(P);
    E=Oscar.Graphs.edges(EdGr);

    graph_bol=true
    for vert in range(1,length(vertices(P)))
        cont_zero=0
        cont=0
        for vert2 in Graphs.neighbors(EdGr,vert)
            cont=cont+1
            if Mukai_form(vert,vert2)==0
                cont_zero=cont_zero+1
            end
        end
        println("For vertex ",vert," we get ",cont_zero,"/",cont, " have scalar product zero")
        if cont!=32
            graph_bol=false
        end
        if cont!=cont_zero
            graph_bol=false
        end
    end

    if graph_bol
        return true;
    end
    
    return false;
end