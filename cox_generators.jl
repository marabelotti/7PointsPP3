using Oscar
using InvertedIndices

include("./helper_functions.jl")

"""

To each section of anticanonical degree one in the Cox ring of the blow up of n=7 points in the projective space of dimension d=3 we can associate a polynomial in the ring CC(t)[x1,x2,x3,x4,x5,x6,x7,y1,y2,y3,y4,y5,y6,y7], as explained in Lemma 7 of the paper with ArXiv number 2208.05258.

In order to reduce the computational effort, in this file we implement functions which returns such polynomials evaluated at x1,x2,x3,x4,x5,x6,x7=1. Knowing the multivariate degree of the sections, we can uniquely reconstruct the original polynomial.
The input of the functions is a 4x7 matrix A with coefficients in CC(t) whose columns represent the projective points in PP3.

The file can be modified to work with different n and d, but not all the sections are obtained in different cases.

"""

n=7;
d=3;
KS, (t,z...)=PolynomialRing(QQ, ["t",["z[$i]" for i in 1:n]...]);

"""
    exceptional_sections()

    
    This function returns a tuple.
    The first output is the 7 element vector of polynomials associated to the exceptional divisors with the x variables evaluated at one. 
    The second output is the vector of multidegrees E_i of these sections.

# Examples
```julia-repl
julia> exceptional_sections();
```

"""


function exceptional_sections(A)
    exc=Vector{fmpq_mpoly}([KS(1) for i in 1:n])
    deg_mon_exc=[[j==i+1 ? 1 : 0 for j in 1:n+1] for i in 1:n]
    return exc,deg_mon_exc
end


"""
    hyperplanes_sections(A::Matrix{fmpq_mpoly})

    This function returns a tuple.
    The first output is the 35 element vector of polynomials associated to hyperplanes with the x variables evaluated at one. 
    The second output is the vector of multidegrees H-E_i-E_j-E_k of these sections.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> hyperplanes_sections(A);
```


"""

function hyperplanes_sections(A::Matrix{fmpq_mpoly})

    @assert size(A)==(d+1,n)

    Hyp=Vector{fmpq_mpoly}([det(matrix(KS,[[z[1:d+1]]; [KS.(A)[:,s] for s in (i,j,k)]])) for i in 1:n for j in i+1:n for k in j+1:n]);
    Hyp=map(y->sub_linear(y,A),Hyp)
    deg_mon_hyp=[[((v==i+1) | (v==j+1) | (v==k+1)) ? 0 : 1 for v in 1:n+1] for i in 1:n for j in i+1:n for k in j+1:n];
    
    return Hyp , deg_mon_hyp
end


"""
    quadrics_sections(A::Matrix{fmpq_mpoly})

    This function returns a tuple.
    The first output is the 42 element vector of polynomials associated to the quadrics with the x variables evaluated at one. 
    The second output is the vector of multidegrees 2H-sum{i!=j} E_i of these sections.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> quadrics_sections(A);
```


"""

function quadrics_sections(A::Matrix{fmpq_mpoly})

    @assert size(A)==(d+1,n)

    Quad=Vector{fmpq_mpoly}()
    for i in 1:n
        for j in 1:n
            if j!=i
                #println("We are doing quadric ",i,j)
                temp=[veronese(z[1:d+1])]
          
                for ind_row in 1:n
                    F = hom(KS, KS, z -> z, [t,A[:,ind_row]...,z[d+2:n]...])
                    if (ind_row!=i) & (ind_row!=j) #i is the excluded point
                        push!(temp,F.(veronese(z[1:d+1])))
                    elseif (ind_row==j)
                        for var in z[1:d+1]
                            push!(temp,F.(derivative.(veronese(z[1:d+1]),var)))
                        end
                    end
                end
                #j is the double point. So we need to look at the derivates to impose multiplicity two
                push!(Quad,det(matrix(KS,temp)))
            end
        end
    end

     
    Quad=map(y->sub_linear(y,A),Quad)

    deg_mon_quad=[[((v==1) | (v==i+1)) ? 2 : (v==j+1 ? 0 : 1) for v in 1:n+1] for i in 1:n for j in ((1:i-1)...,(i+1:n)...)]
    
        return Quad , deg_mon_quad
 end


"""
    cubics_sections(A::Matrix{fmpq_mpoly})

    This function returns a tuple.
    The first output is the 35 element vector of polynomials associated to the cubics with the x variables evaluated at one. 
    The second output is the vector of multidegrees 3H-E_i-E_j-E_k-2 sum{s!=i,j,k} E_s of these sections.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> cubics_sections(A);
```


"""

function cubics_sections(A::Matrix{fmpq_mpoly})

    @assert size(A)==(d+1,n)

    Cubic=Vector{fmpq_mpoly}()
    for i in 1:7
        for j in 1:i-1
            for k in 1:j-1

                change_matrix=matrix(KS,A[:,Not(i,j,k)])
                points_not_double=matrix(KS,A[:,[i,j,k]])
               
                B=cofactor(change_matrix)*points_not_double

                const_mat=zero_matrix(KS,d,d+1)
                const_mat[1,1]=B[1,3]*B[2,3]*B[3,3]
                const_mat[2,1]=B[1,2]*B[2,2]*B[3,2]
                const_mat[3,1]=B[1,1]*B[2,1]*B[3,1]

                const_mat[1,2]=B[1,3]*B[2,3]*B[4,3]
                const_mat[2,2]=B[1,2]*B[2,2]*B[4,2]
                const_mat[3,2]=B[1,1]*B[2,1]*B[4,1]

                const_mat[1,3]=B[1,3]*B[3,3]*B[4,3]
                const_mat[2,3]=B[1,2]*B[3,2]*B[4,2]
                const_mat[3,3]=B[1,1]*B[3,1]*B[4,1]
            
                const_mat[1,4]=B[2,3]*B[3,3]*B[4,3]
                const_mat[2,4]=B[2,2]*B[3,2]*B[4,2]
                const_mat[3,4]=B[2,1]*B[3,1]*B[4,1]

                const_vec=[(-1)^(i-1)*det(const_mat[:, [(1:i-1)...,(i+1:d+1)...]]) for i in 1:d+1]
                cubic_temp=sum([const_vec[i]*prod(z[Not(d+2-i,d+2:n...)]) for i in 1:d+1])
                
                temp = cofactor(change_matrix)*z[1:d+1]
                F = hom(KS, KS, z -> z, [t,temp[1:d+1]...,z[d+2:n]...])
                
                push!(Cubic,F(cubic_temp))
            end
        end
    end

     
    Cubic=map(y->sub_linear(y,A),Cubic)

    deg_mon_cubic=[[(v==1 ? 3 : ((v==i+1 || v==j+1 || v==k+1) ? 2 : 1)) for v in 1:n+1] for i in 1:n for j in 1:i-1 for k in 1:j-1]
    
        return Cubic , deg_mon_cubic
end

"""
    quartic_sections(A::Matrix{fmpq_mpoly},quad::Vector{fmpq_mpoly},stran::Vector{fmpq_mpoly})

    The input of this function are the usual matrix A and the first outputs of quadric_sections and strange_sections.
    This function returns a tuple.
    The first output is the 7 element vector of polynomials associated to the quartics with the x variables evaluated at one. 
    The second output is the vector of multidegrees 4H-3E_i-sum{i!=j} E_j of these sections.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> quartics_sections(A);
```


"""

function quartics_sections(A::Matrix{fmpq_mpoly},quad::Vector{fmpq_mpoly},stran::Vector{fmpq_mpoly})

    @assert size(A)==(d+1,n)

    Quart=Vector{fmpq_mpoly}()
    
    label=[(i,j) for i in 1:n for j in (1:i-1...,i+1:n...)]
    
    ind_0=[2,1,1,1,1,1,1]
    ind_j=[3,3,2,3,4,5,6]
    
    for i in 1:7
        println("We are doing quartic ",i)
        flush(stdout)
        Q0i=quad[findfirst(x->x==(ind_0[i],i), label)];
        Q0j=quad[findfirst(x->x==(ind_0[i],ind_j[i]), label)];
        Qji=quad[findfirst(x->x==(ind_j[i],i), label)];
        S0=stran[1];
        S1=stran[2];
        S2=stran[3];
        eq_mono=[S0*Q0i,S1*Q0i,S2*Q0i,Q0j*Qji];
        list_mon=[monomials(S0*Q0i)...,monomials(S1*Q0i)...,monomials(S2*Q0i)...,monomials(Q0j*Qji)...];
        mono=Set([mon for mon in list_mon]);

        b=[]
        for cont in range(1,7)
            if cont!=ind_0[i] && cont!=i 
                    push!(b,[(v==ind_0[i]) ? 3 : (v==cont ? 1 : 0) for v in 1:n])
            end
        end
        
        M_tar=zero_matrix(FractionField(KS),4,6);

        for (el,col) in zip(b,range(1,6))
            for (eq,row) in zip(eq_mono,range(1,4))
                M_tar[row,col]=coeff(eq,[2,3,4,5,6,7,8],el)
            end
        end

        M_tar=transpose(M_tar);
        ker_M_tar=kernel(M_tar)[2];

        push!(Quart,numerator(sum([eq_mono[i]*ker_M_tar[i,1] for i in 1:4])));
    end

    deg_mon_quartic=[[(v==1 ? 4 : ((v==i+1) ? 1 : 2)) for v in 1:n+1] for i in 1:n]
    
    return Quart , deg_mon_quartic
end

"""
    strange_sections(A::Matrix{fmpq_mpoly},A_plus)

    
    This function returns a tuple.
    The first output is the 3 element vector of polynomials associated to the sections of -K/2 with the x variables evaluated at one.
    The second output is the vector of multidegrees 2H-E_1-E_2-E_3-E_4-E_5-E_6-E_7 of these sections.

# Examples
```julia-repl
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia>A_plus=Matrix{fmpq_mpoly}([KS(0) t^2 KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) t^3; KS(0) KS(0) KS(0) KS(0) KS(0) t^4 KS(0) KS(0) KS(0) KS(0); KS(1) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(1) KS(0) KS(0) KS(0); KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(0) KS(1) KS(0)])
julia> strange_sections(A,A_plus);
```


"""

function strange_sections(A::Matrix{fmpq_mpoly},A_plus::Matrix{fmpq_mpoly})

    @assert size(A)==(d+1,n)

    Strange=Vector{fmpq_mpoly}()
    
    mat_lin_ind=zero_matrix(KS,10,3)
    
    for i in 1:3
        temp=[veronese(z[1:d+1])]
        for ind_row in 1:n
            push!(temp,veronese(A[:,ind_row]))
        end
        
        for ind_row in 1:2
            push!(temp,A_plus[2*i-2+ind_row,:])
        end

        push!(Strange,det(matrix(KS,temp)))
        
        for ind in 1:10
            mat_lin_ind[ind,i]=coeff(det(matrix(KS,temp)),veronese(z[1:d+1])[ind])
        end
    end
     
    Strange=map(y->sub_linear(y,A),Strange)

    deg_mon_stran=[[2 1 1 1 1 1 1 1] for i in 1:3]
    
    return Strange , deg_mon_stran, (rank(mat_lin_ind)==3)
end
