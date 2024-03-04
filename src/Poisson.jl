function Coulomb(vert_i::Vertex{T}, vert_j::Vertex{T}) where{T}
    return one(T) / norm(vert_i.p - vert_j.p)
end

function ∂ₙᵢCoulomb(vert_i::Vertex{T}, vert_j::Vertex{T}) where{T}
    return - dot(vert_i.p - vert_j.p, vert_i.n) / norm(vert_i.p - vert_j.p)^3
end

function Poission_LHS()

end

function Poission_RHS()

end