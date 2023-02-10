using Oscar

"""
    leading_t(eq::fmpq_mpoly)

            The input eq must be a polynomial in some multivariate polynomial ring R::FmpqMPolyRing.
It returns the coefficent of smallest degree of eq when seen as a polynomial in one variable, specifically the first variable appearing in the definition of R.

# Examples
```julia-repl
julia> KS, (t,z0,z1)=PolynomialRing(QQ, ["t","z0","z1"]);
julia> eq=t^2*z0*z1+t*z0+t*z1;
julia> leading_t(eq)
z0 + z1
```

"""
function leading_t(eq::fmpq_mpoly)

    if(length(eq)==0)
        return eq,true;
    end
    
    sort_terms!(eq)
    deg_lead=degree(term(eq,length(eq)),t)
    lead=coeff(eq,[1],[deg_lead])
    monomial=ismonomial(lead)
    return lead;
end


"""

    leading_t_vector(eq_vec::Vector{fmpq_mpoly})

    The input is a vector eq_vect of polynomials in some multivariate polynomial ring R::FmpqMPolyRing. It returns a tuple.
The first element of the tuple is a vector of polynomials obtained by applying leading_t to every element of eq_vec.
The second element is a boolean which is true if the leading coefficient of every element is monomial and false if there exists one for which it is not.


"""

function leading_t_vect(eq_vect::Vector{fmpq_mpoly})
    return  (leading_t.(eq_vect)), all((length.(leading_t.(eq_vect))).<2);
end


"""

    veronese(v::Vector{Any})

    This encodes the veronese map on the vector v.

"""


function veronese(v::Vector{fmpq_mpoly})
    return [v[i]*v[j] for i in 1:length(v) for j in i:length(v)]
end


"""

    sub_linear(A::Matrix{fmpq_mpoly},eq::fmpq_mpoly)

    The input of this function are a matrix A of dimension n1xn2 and a polynomial eq::fmpq_mpoly in a multivariate polynomial ring R.
    The output is a polynomial of R obtained by substituting the variables x_{2},..,x_{n1+1} with the columns of A*[x_{2},..,x{n1+1}] 

# Examples
```julia-repl
julia>KS, (t,z0,z1,z2,z3,z4)=PolynomialRing(QQ, ["t","z0","z1","z2","z3","z4"]);
julia>eq=t*z0+t^2*z4;
julis> A=[0 2;1 2; 3 4];
julia> sublinear(A,eq)
t^2*z4 + 2*t*z1
```

"""

function sub_linear(eq::fmpq_mpoly,A::Matrix{fmpq_mpoly})
    linear_forms=Vector{fmpq_mpoly}(A*gens(eq.parent)[2:1+size(A)[2]])
    sub = hom(eq.parent, eq.parent, z -> z, [eq.parent(t),linear_forms..., gens(eq.parent)[size(A)[1]+2:end]...])
    return sub(eq)
end


"""

    cofactor(A::AbstractAlgebra.Generic.MatSpaceElem{fmpq_mpoly})

    The input of this function is a matrix A.
    The output is the adjugate matrix of A. 


"""


function cofactor(A::AbstractAlgebra.Generic.MatSpaceElem{fmpq_mpoly})
   ax = axes(A)
   out = zero_matrix(KS,size(A)...)
   for col in ax[1]
       for row in ax[2]
            out[row, col] = (-1)^(col + row) *det(A[[(1:col-1)...,(col+1:size(A,1))...],[(1:row-1)...,(row+1:size(A,2))...]])
       end
   end
   return out
end


"""

    zero_check(eq::Vector{fmpq_mpoly})

    The input of this function is a vectors of polynomials.
    The output is true if every element in the Vector is different from the zero polynomial, false otherwise.


"""

function zero_check(eq)
    return all(eq .!=KS(0))
end

"""

    strange_diff(eq::Vector{fmpq_mpoly})

    The input of this function is a vector of size three.
    The output is true if there are two elements in the vector who share the same leading monomial in the x, false otherwise.


"""


function strange_diff(stra_mon)

    @assert length(stra_mon)==3;
    
    bol=true
    mon=[numerator(monomial(el,1)//t^(degree(el,t))) for el in stra_mon]

    for i in 1:2
        for j in i+1:3
            if (mon[i]==mon[j]) 
                println("Two monomials are the same")
                bol=false
            end
        end
    end
    
    return bol;   
end


"""

    degree_quadratic(deg_mon::Vector{Vector{Int64}, equations::Vector{fmpq_mpoly})

    This function has two functions.
    The first output is a vector of couples producing those multidegrees and their indices.
    The second output is a vector of all possible multidegrees of anticanonical degree two. 
"""

function degree_quadratic(deg_mon,equations)
    couples=[[equations[i],equations[j],i,j] for i in 1:129 for j in i:129]
    deg_couples=[deg_mon[i].+deg_mon[j] for i in 1:129 for j in i:129]
    return couples,deg_couples
end
