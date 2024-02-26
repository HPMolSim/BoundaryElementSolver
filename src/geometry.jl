# the operations about triangular meshes

@inline function tri_mean(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    return (a + b + c) ./ 3, 0.5 * norm(normal)
end

@inline function tri_area(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    normal = cross(b - a, c - a)
    return 0.5 * norm(normal)
end