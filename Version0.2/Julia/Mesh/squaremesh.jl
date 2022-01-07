
function TriSquareMesh(
    m::Int,
    n::Int = m,
)::Tuple{Array{Float64,2},Array{Int,2}}

    p = [[x; y] for x = 0.0:1.0/Float64(m):1.0, y = 0.0:1.0/Float64(n):1.0][:]
    t = [
        [
            [i + (j - 1) * (m + 1) i + (j - 1) * (m + 1) + 1 i + j * (m + 1)],
            [i + (j - 1) * (m + 1) + 1 i + j * (m + 1) + 1 i + j * (m + 1)],
        ] for i = 1:m, j = 1:n
    ]

    return hcat(p...), hcat(hcat(t'...)...)
end

function QuadSquareMesh(
    m::Int,
    n::Int = m,
)::Tuple{Array{Float64,2},Array{Int,2}}

    p = [[x; y] for x = 0.0:1.0/Float64(m):1.0, y = 0.0:1.0/Float64(n):1.0][:]
    t = [ [i + (j - 1) * (m + 1)    i + (j - 1) * (m + 1) + 1  i + j * (m + 1) + 1   i + j * (m + 1)]
          for j = 1:m, i = 1:n ]

    return hcat(p...), hcat(t'...)
end

function SquareMesh(
    m::Int,
    n::Int,
    elemtype::Int,
)::Tuple{Array{Float64,2},Array{Int,2}}

    if (elemtype==0)
        p, t = TriSquareMesh(m,n)
    else
        p, t = QuadSquareMesh(m,n)
    end

    return p, t
end
