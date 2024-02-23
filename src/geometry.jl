# the operations about triangular meshes

@inline function triangle(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    normal = cross(b - a, c - a)
    return (a + b + c) ./ 3, normalize(normal), 0.5 * norm(normal)
end