using Oscar
using Combinatorics
                                  
S_t, (t,z0,z1,z2,z3,z4,z5,z6,e0,e1,e2,e3,e4,e5,e6,H012,H013,H014,H015,H016,H023,H024,H025,H026,H034,H035,H036,H045,H046,H056,H123,H124,H125,H126,H134,H135,H136,H145,H146,H156,H234,H235,H236,H245,H246,H256,H345, H346,H356,H456,Q01,Q02,Q03,Q04,Q05,Q06,Q10,Q12,Q13,Q14,Q15,Q16,Q20,Q21,Q23,Q24,Q25,Q26,Q30,Q31,Q32,Q34,Q35,Q36,Q40,Q41,Q42,Q43,Q45,Q46,Q50,Q51,Q52,Q53,Q54,Q56,Q60,Q61,Q62,Q63,Q64,Q65,C210,C310,C320,C321,C410,C420, C421,C430,C431,C432,C510,C520,C521,C530,C531,C532,C540,C541,C542,C543,C610,C620,C621,C630,C631,C632,C640,C641,C642,C643,C650,C651,C652,C653,C654,Qr0,Qr1,Qr2,Qr3,Qr4,Qr5,Qr6,S0,S1,S2)= PolynomialRing(QQ, ["t","z0","z1","z2","z3","z4","z5","z6","e0","e1","e2","e3","e4","e5","e6","H012","H013","H014","H015","H016","H023","H024","H025","H026","H034","H035","H036","H045","H046","H056","H123","H124","H125","H126","H134","H135","H136","H145","H146","H156","H234","H235","H236","H245","H246","H256","H345","H346","H356","H456","Q01","Q02","Q03","Q04","Q05","Q06","Q10","Q12","Q13","Q14","Q15","Q16","Q20","Q21","Q23","Q24","Q25","Q26","Q30","Q31","Q32","Q34","Q35","Q36","Q40","Q41","Q42","Q43","Q45","Q46","Q50","Q51","Q52","Q53","Q54","Q56","Q60","Q61","Q62","Q63","Q64","Q65","C210","C310","C320","C321","C410","C420","C421","C430","C431","C432","C510","C520","C521","C530","C531","C532","C540","C541","C542","C543","C610","C620","C621","C630","C631","C632","C640","C641","C642","C643","C650","C651","C652","C653","C654","Qr0","Qr1","Qr2","Qr3","Qr4","Qr5","Qr6","S0","S1","S2"]);

include("./cox_generators.jl")

A_plus_default=Matrix{fmpq_mpoly}([KS(0) t^2 KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) t^3; KS(0) KS(0) KS(0) KS(0) KS(0) t^4 KS(0) KS(0) KS(0) KS(0); KS(1) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(1) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(1) KS(0)])

"""
    generalposition_monomeric_check(A::Matrix{fmpq_mpoly},A_plus::Matrix{fmpq_mpoly}=A_plus_default)

    The input of this function is just a matrix A with entries in KS=CC[t], and a matrix A_plus necessary to compute the strange_sections to which we gave a default value.
    The output are:
    1) A boolean  which is true if and only if the seven points given by the columns of P are in general position as described in Section 3 of the article and the initial forms are monomeric
    2)If bool is true, a vector of integers deg_mon which gives the multidegree of the 129 generators
    3)If bool is true, a vector of monomials in KS which are the initial forms of the images of the generators through phi(t) but with the x coordinates evaluated at one and the y called z. This is just a raw version of the initial forms which we will obtained later by calling the function initial_forms_of_sections


# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> bool,deg_mon,equations_mon=generalposition_monomeric_check(A);
```

"""

function generalposition_monomeric_check(A::Matrix{fmpq_mpoly},A_plus::Matrix{fmpq_mpoly}=A_plus_default)
    
    exc, deg_mon_exc=exceptional_sections(A);
    exc_mon,bol_exc=leading_t_vect(exc);

    hyps, deg_mon_hyp= hyperplanes_sections(A);
    if !zero_check(hyps)
        println("One of the hyps is zero")
        return false;
    end
    hyps_mon,bol_hyps= leading_t_vect(hyps);
    if !bol_hyps
        println("The hyperplanes are not monomeric")
        return false;
    end
    println("Hyperplanes check ")
    flush(stdout)

    quad,deg_mon_quad = quadrics_sections(A);
    if !zero_check(quad)
        println("One of the quad is zero")
        return false;
    end
    quad_mon,bol_quad= leading_t_vect(quad);
    if !bol_quad
        println("The quadrics are not monomeric")
        return false;
    end
    println("Quad check ")
    flush(stdout)

    cub,deg_mon_cub = cubics_sections(A);
    if !zero_check(cub)
        println("One of the cub is zero")
        return false;
    end
    cub_mon,bol_cub= leading_t_vect(cub);
    if !bol_cub
        println("The cubics are not monomeric")
        return false;
    end
    println("Cub check ")
    flush(stdout)

    stra,deg_mon_strange,bol_strange=strange_sections(A,A_plus);
    if !zero_check(stra)
        println("One of the strange is zero")
        return false;
    end
    stra_mon,bol_strange=leading_t_vect(stra);
    if  !bol_strange
        println("The strange sections are not monomeric")
        return false;
    end
    if !strange_diff(stra_mon)
        println("Two initial forms of the strange sections are the same. Change A_bol!")
        return false;
    end
    println("Strange check")
    flush(stdout)
     
    quart,deg_mon_quart = quartics_sections(A,quad,stra);
    if !zero_check(quart)
        println("One of the quartics is zero")
        return false;
    end
    quart_mon,bol_quart= leading_t_vect(quart);
    if !bol_quart
        println("The quartics are not monomeric")
        return false;
    end
    println("Quartics check")
    flush(stdout)

    deg_mon=[deg_mon_exc...,deg_mon_hyp...,deg_mon_quad...,deg_mon_cub...,deg_mon_quart...,deg_mon_strange...]
    equations_mon=[exc_mon...,hyps_mon...,quad_mon...,cub_mon...,quart_mon..., stra_mon...];

    return true,deg_mon,equations_mon;
end


"""
    degree_two_check(deg_mon::Vector{Vector{Int64}},equations_mon::Vector{fmpq_mpoly})

    The inputs of this function are a vector deg_mon of the multidegrees of the 129 generators of the Cox Ring of the blow up of 7 points in PP3 and a vector of 129 monomials in CC[z1,z2,z3,z4,z5,z6,z7]  which are the initial forms of the generators but with the x variables evaluated at 1 and the y are called z.  
    The output is a boolean which is true if and only if the Hilbert function of the algebra generated by the initial forms of the images of phi(t) is equal to the Hilbert function of the the Cox ring in degree two.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> bool,deg_mon,equations_mon=generalposition_monomeric_check(A);
julia> degree_two_check(deg_mon,equations_mon)
true
```

"""

function degree_two_check(deg_mon::Vector{Vector{Int64}},equations_mon::Vector{fmpq_mpoly})
    couples,deg_couples=degree_quadratic(deg_mon,equations_mon);

    vars_S=gens(S_t)[n+2:end];
    same_ring=hom(KS,S_t,z->z,[gens(S_t)[1:8]...])

    dimensions=[2,4,2,2,4,7,4,2,2,2,4,2,4,2,2]
    orbit_degreetwo_polymake=[[1, 1, 0, 0, 1, 1, 1, 1], [2, 1, 2, 1, 1, 1, 1, 1], [2, 2, 2, 0,1, 1, 1, 1], [3, 3, 2, 2, 1, 1, 1, 1], [3, 2, 2, 2, 2, 1, 1, 1], [4, 2, 2, 2, 2, 2, 2, 2], [4, 3, 1, 2, 2, 2, 2, 2], [3, 2, 0, 1, 2, 2, 2, 2], [4, 2, 1, 1, 3, 3, 2, 2], [5, 1, 2, 2, 3, 3, 3, 3], [5, 2, 2, 2, 2, 3, 3, 3], [5, 4, 2, 3, 2, 2, 2, 2], [6, 3, 2, 3, 3, 3, 3, 3], [6, 4, 2, 2, 3, 3, 3, 3], [7, 3, 4, 4, 3, 3, 3, 3]];

    bol_hilb=true
    for (deg,cont_deg) in zip(orbit_degreetwo_polymake,range(1,length(orbit_degreetwo_polymake)))
        println("--------------------------- ",cont_deg," out of ",length(orbit_degreetwo_polymake),". Check of degree ", deg)
        #println("There are ",length(multiset_permutations(deg[2:8],7)) ," permutations")
        for deg_perm in multiset_permutations(deg[2:8],7)
            pushfirst!(deg_perm,deg[1])
            monomials=[]
            
            for (cou,i) in zip(couples,range(1,length(couples)))
                if deg_couples[i]==deg_perm
                    push!(monomials,same_ring(cou[2])*same_ring(cou[1]));
                end
            end
            if length(monomials)!=1
                monomials_mon=Set()
                for el in monomials
                    push!(monomials_mon,monomial(el,1))
                end
                if length(monomials_mon)!=dimensions[cont_deg]
                    println(deg_perm)
                    println(monomials_mon)
                    bol_hilb=false
                end
            end
        end
    end
    
    if !bol_hilb
        println("Not Khovanskii cause Hilbert function of initial algebra doesn't match")
        return false;
    end
    #Now we know it is Khovanskii. We now compute the change of variables z_i -> y_i/x_i and the homogenization.     

    return true
end

"""
    initial_forms_of_sections(equations_mon::Vector{fmpq_mpoly})

    The input of this function is a vector of 129 monomials in CC[z1,z2,z3,z4,z5,z6,z7]  which are the initial forms of the generators but with the x variables evaluated at 1 and the y are called z.  
    The output is the vector of the 129 initial forms as polynomials in CC[x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6].

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> bool,deg_mon,equations_mon=generalposition_monomeric_check(A);
julia> initial_form=degree_two_check(equations_mon);
```

"""

function initial_forms_of_sections(equations_mon::Vector{fmpq_mpoly})
    exc_mon=equations_mon[1:7]
    hyp_mon=equations_mon[8:7+35]
    quad_mon=equations_mon[43:42+42]
    cub_mon=equations_mon[85:84+35]
    quart_mon=equations_mon[120:119+7]
    stra_mon=equations_mon[127:end]
    
    R, (t,x1,x2,x3,x4,x5,x6,x7,y1,y2,y3,y4,y5,y6,y7) = PolynomialRing(QQ, ["t","x1","x2","x3","x4","x5","x6","x7","y1","y2","y3","y4","y5","y6","y7"])
    x = gens(R)[2:8];
    y = gens(R)[9:end];

    exc_mon_R=[x[i] for i in  1:7];
    hyps_mon_R,quad_mon_R,cub_mon_R,quart_mon_R,stra_mon_R=[],[],[],[],[];

    Toxy=hom(KS,FractionField(R),z->z,[t,y1//x1,y2//x2,y3//x3,y4//x4,y5//x5,y6//x6,y7//x7])

    cont=1;
    for i in 1:7
        for j in i+1:7
            for k in j+1:7
                push!(hyps_mon_R,numerator(Toxy(hyps_mon[cont])*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]//(x[i]*x[j]*x[k])));
                cont=cont+1;
            end
        end
    end

    cont=1;
    for i in 1:7
        for j in 1:7
            if j!=i
                push!(quad_mon_R,numerator(Toxy(quad_mon[cont])*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*x[i]//x[j]));
                cont=cont+1;
            end
        end
    end

    cont=1;
    for i in 1:7
        for j in 1:i-1
            for k in 1:j-1
                push!(cub_mon_R,numerator(Toxy(cub_mon[cont])*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*x[i]*x[j]*x[k]));
                cont=cont+1;
            end
        end
    end

    for i in 1:7
        push!(quart_mon_R,numerator(Toxy(quart_mon[i])*x[1]^2*x[2]^2*x[3]^2*x[4]^2*x[5]^2*x[6]^2*x[7]^2//x[i]));
    end

    for i in 1:3
        push!(stra_mon_R,numerator(Toxy(stra_mon[i])*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]));
    end

    return exc_mon_R,hyps_mon_R,quad_mon_R,cub_mon_R,quart_mon_R,stra_mon_R;
end

"""
    matrix_of_exponents(initial_forms::Vector{fmpq_mpoly})
    
    The input of this function is a vector of the 129 initial forms as polynomials in CC[x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6].
    The output is just a 129x14 matrix of the exponents of these initial forms.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> bool,deg_mon,equations_mon=generalposition_monomeric_check(A);
julia> initial_forms=degree_two_check(equations_mon);
julia>matrix_of_exponents(initial_forms);
```

"""

function matrix_of_exponents(initial_forms::Vector{fmpq_mpoly})
    Ring=equations_mon[1].parent;
    M=zero_matrix(ZZ,129,14);
    for (eq,row) in zip(initial_forms,range(1,129))
        for (var,col) in zip(gens(Ring)[2:end],range(1,14))
            M[row,col]=Oscar.degree(eq,var)
        end
    end
    return M; 
end

