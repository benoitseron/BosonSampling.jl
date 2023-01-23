"""
    gram_matrix_parametrized(α,β,γ,ϕ)

Using minimal elements to parametrize a n = 3 gram matrix;
we look at appendix 1 of 1902.11039 to use the three mode parametrization of the states, and we then generate a gram matrix from this.
"""
function gram_matrix_parametrized(α,β,γ,ϕ)

    ab = cos(β)
    ac = cos(γ)
    bc = ab*ac + exp(1im * ϕ) *sin(β)*sin(γ)*sin(α)

    G = [1 ab ac; conj(ab) 1 bc; conj(ac) conj(bc) 1]

end
