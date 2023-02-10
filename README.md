# 7PointsPP3
This repository contains the auxiliary code to the paper "DISCRETE GEOMETRY OF COX RINGS OF BLOW-UPS OF $\mathbb{P}^3$", written by Mara Belotti and Marta Panizzut.

### About the paper
The Arxiv identifier 2208.05258.  We prove quadratic generation for the ideal of the Cox ring of
the blow-up of $\mathbb{P}^3$ at 7 points, solving a conjecture of Lesieutre and Park. To
do this we compute Khovanskii bases, implementing techniques which proved
successful in the case of Del Pezzo surfaces. Such bases give us degenera-
tions to toric varieties whose associated polytopes encode toric degenerations
with respect to all projective embeddings. We study the edge-graphs of these
polytopes and we introduce the Mukai edge graph.


### About the code
The helper_functions.jl and cox_generators.jl contains auxiliary code for the main functions in is_khovanskii.jl.
This is a quick example on how to use the code in a julia terminal:

```julia
julia> KS, (t,z...)=PolynomialRing(QQ, ["t",["z[$i]" for i in 1:7]...]);
julia> A=[t^59 t^44 t^79 t^20 t^12 t^81 t^36; t^8 t^72 t^49 t^39 t^58 t^23 t^64; t^44 t^58 t^12 t^52 t^57 t^49 t^51; t^25 t^23 t^60 t^72 t^45 t^51 t^6];
julia> bool,_=is_khovanskii(A);
julia> bool
true

```
The necessary code for the proof of Theorem 12 is contained in the tutorial.jl jupyter notebook. The file Lemma22.ipynb contains the code necessary for the proof of Corollary 22.
