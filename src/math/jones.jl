# Copyright (c) 2015-2017 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

abstract type AbstractJonesMatrix <: AbstractMatrix{ComplexF64} end

@doc doc"""
    JonesMatrix

This type represents a 2x2 complex Jones matrix.

```math
    \begin{pmatrix}
        j_{xx} & j_{xy} \\
        j_{yx} & j_{yy} \\
    \end{pmatrix}
```
"""
struct JonesMatrix <: AbstractJonesMatrix
    xx::ComplexF64
    xy::ComplexF64
    yx::ComplexF64
    yy::ComplexF64
end

@doc doc"""
    DiagonalJonesMatrix

This type represents a Jones matrix that is diagonal.

```math
    \begin{pmatrix}
        j_{xx} & 0 \\
        0 & j_{yy} \\
    \end{pmatrix}
```

These matrices are used to represent the complex gains of each antenna
without accounting for the off-diagonal polarization leakage terms.
"""
struct DiagonalJonesMatrix <: AbstractJonesMatrix
    xx::ComplexF64
    yy::ComplexF64
end

@doc doc"""
    HermitianJonesMatrix

This type represents a Jones matrix that is Hermitian.

```math
    \begin{pmatrix}
        j_{xx} & j_{xy} \\
        j_{xy}^* & j_{yy} \\
    \end{pmatrix}
```

These matrices are useful for representing the $xx$, $xy$, $yx$, and $yy$ flux,
because $xx$ and $yy$ are constrained to be real while $xy$ and $yx$ are complex
conjugates.
"""
struct HermitianJonesMatrix <: AbstractJonesMatrix
    xx::Float64
    xy::ComplexF64
    yy::Float64
end

# Abstract array interface

Base.size(::AbstractJonesMatrix) = (2, 2)

@noinline function checkbounds(J, idx)
    if idx < 1 || idx > 4
        throw(BoundsError(J, idx))
    end
end

function Base.getindex(J::JonesMatrix, idx::Int)
    @boundscheck checkbounds(J, idx)
    if idx == 1
        J.xx
    elseif idx == 2
        J.yx
    elseif idx == 3
        J.xy
    else
        J.yy
    end
end

function Base.getindex(J::DiagonalJonesMatrix, idx::Int)
    @boundscheck checkbounds(J, idx)
    if idx == 1
        J.xx
    elseif idx == 4
        J.yy
    else
        zero(ComplexF64)
    end
end

function Base.getindex(J::HermitianJonesMatrix, idx::Int)
    @boundscheck checkbounds(J, idx)
    if idx == 1
        J.xx
    elseif idx == 2
        conj(J.xy)
    elseif idx == 3
        J.xy
    else
        J.yy
    end
end

Base.IndexStyle(::Type{<:AbstractJonesMatrix}) = IndexLinear()

Base.zero(::Type{JonesMatrix}) = JonesMatrix(0, 0, 0, 0)
Base.one( ::Type{JonesMatrix}) = JonesMatrix(1, 0, 0, 1) # the identity matrix
Base.rand(::Type{JonesMatrix}) = JonesMatrix(rand(ComplexF64), rand(ComplexF64),
                                             rand(ComplexF64), rand(ComplexF64))

Base.zero(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(0, 0)
Base.one( ::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(1, 1) # the identity matrix
Base.rand(::Type{DiagonalJonesMatrix}) = DiagonalJonesMatrix(rand(Float64), rand(Float64))

Base.zero(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(0, 0, 0)
Base.one( ::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(1, 0, 1) # the identity matrix
Base.rand(::Type{HermitianJonesMatrix}) = HermitianJonesMatrix(rand(Float64), rand(ComplexF64),
                                                               rand(Float64))

#function JonesMatrix(mat::Matrix)
#    size(mat) == (2,2) || throw(DimensionMismatch("A Jones matrix must be 2x2."))
#    JonesMatrix(mat[1,1], mat[1,2], mat[2,1], mat[2,2])
#end

for op in (:+, :-)
    @eval function Base.$op(J1::JonesMatrix, J2::JonesMatrix)
        JonesMatrix($op(J1.xx, J2.xx),
                    $op(J1.xy, J2.xy),
                    $op(J1.yx, J2.yx),
                    $op(J1.yy, J2.yy))
    end
    @eval function Base.$op(J1::DiagonalJonesMatrix, J2::DiagonalJonesMatrix)
        DiagonalJonesMatrix($op(J1.xx, J2.xx),
                            $op(J1.yy, J2.yy))
    end
    @eval function Base.$op(J1::HermitianJonesMatrix, J2::HermitianJonesMatrix)
        HermitianJonesMatrix($op(J1.xx, J2.xx),
                             $op(J1.xy, J2.xy),
                             $op(J1.yy, J2.yy))
    end
end

for op in (:*, :/)
    @eval function Base.$op(J::JonesMatrix, a::Number)
        JonesMatrix($op(J.xx, a), $op(J.xy, a), $op(J.yx, a), $op(J.yy, a))
    end
    @eval function Base.$op(J::DiagonalJonesMatrix, a::Number)
        DiagonalJonesMatrix($op(J.xx, a), $op(J.yy, a))
    end
    @eval function Base.$op(J::HermitianJonesMatrix, a::Real)
        HermitianJonesMatrix($op(J.xx, a), $op(J.xy, a), $op(J.yy, a))
    end
    @eval function Base.$op(J::HermitianJonesMatrix, a::Number)
        JonesMatrix($op(J.xx, a), $op(J.xy, a), $op(conj(J.xy), a), $op(J.yy, a))
    end
end

Base.:*(a::Number, J::AbstractJonesMatrix) = J*a
Base.:/(a::Number, J::AbstractJonesMatrix) = a*inv(J)
Base.:\(a::Number, J::AbstractJonesMatrix) = J/a

function Base.:*(J1::JonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.yx, J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.yx, J1.yx*J2.xy + J1.yy*J2.yy)
end

function Base.:*(J1::JonesMatrix, J2::DiagonalJonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xy*J2.yy, J1.yx*J2.xx,J1.yy*J2.yy)
end

function Base.:*(J1::DiagonalJonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx*J2.xx,J1.xx*J2.xy, J1.yy*J2.yx,J1.yy*J2.yy)
end

function Base.:*(J1::DiagonalJonesMatrix, J2::DiagonalJonesMatrix)
    DiagonalJonesMatrix(J1.xx*J2.xx, J1.yy*J2.yy)
end

function Base.:*(J1::HermitianJonesMatrix, J2::HermitianJonesMatrix)
    JonesMatrix(J1.xx *J2.xx + J1.xy*J2.xy', J1.xx *J2.xy + J1.xy*J2.yy,
                J1.xy'*J2.xx + J1.yy*J2.xy', J1.xy'*J2.xy + J1.yy*J2.yy)
end

function Base.:*(J1::HermitianJonesMatrix, J2::JonesMatrix)
    JonesMatrix(J1.xx *J2.xx + J1.xy*J2.yx, J1.xx *J2.xy + J1.xy*J2.yy,
                J1.xy'*J2.xx + J1.yy*J2.yx, J1.xy'*J2.xy + J1.yy*J2.yy)
end

function Base.:*(J1::JonesMatrix, J2::HermitianJonesMatrix)
    JonesMatrix(J1.xx*J2.xx + J1.xy*J2.xy', J1.xx*J2.xy + J1.xy*J2.yy,
                J1.yx*J2.xx + J1.yy*J2.xy', J1.yx*J2.xy + J1.yy*J2.yy)
end

Base.:\(J1::AbstractJonesMatrix, J2::AbstractJonesMatrix) = inv(J1)*J2
Base.:/(J1::AbstractJonesMatrix, J2::AbstractJonesMatrix) = J1*inv(J2)

Base.conj(J::JonesMatrix) = JonesMatrix(conj(J.xx), conj(J.xy), conj(J.yx), conj(J.yy))
Base.conj(J::DiagonalJonesMatrix) = DiagonalJonesMatrix(conj(J.xx), conj(J.yy))
Base.conj(J::HermitianJonesMatrix) = HermitianJonesMatrix(J.xx, conj(J.xy), J.yy)

LinearAlgebra.transpose(J::JonesMatrix) = JonesMatrix(J.xx, J.yx, J.xy, J.yy)
LinearAlgebra.transpose(J::DiagonalJonesMatrix) = J
LinearAlgebra.transpose(J::HermitianJonesMatrix) = conj(J)

LinearAlgebra.ctranspose(J::JonesMatrix) = JonesMatrix(conj(J.xx), conj(J.yx), conj(J.xy), conj(J.yy))
LinearAlgebra.ctranspose(J::DiagonalJonesMatrix) = conj(J)
LinearAlgebra.ctranspose(J::HermitianJonesMatrix) = J

LinearAlgebra.adjoint(J::JonesMatrix) = JonesMatrix(conj(J.xx), conj(J.yx), conj(J.xy), conj(J.yy))
LinearAlgebra.adjoint(J::DiagonalJonesMatrix) = conj(J)
LinearAlgebra.adjoint(J::HermitianJonesMatrix) = J

LinearAlgebra.det(J::JonesMatrix) = J.xx*J.yy - J.xy*J.yx
LinearAlgebra.det(J::DiagonalJonesMatrix) = J.xx*J.yy
LinearAlgebra.det(J::HermitianJonesMatrix) = J.xx*J.yy - J.xy*conj(J.xy)

function LinearAlgebra.inv(J::JonesMatrix)
    1/det(J) * JonesMatrix(J.yy, -J.xy, -J.yx, J.xx)
end

function LinearAlgebra.inv(J::DiagonalJonesMatrix)
    DiagonalJonesMatrix(1/J.xx, 1/J.yy)
end

function LinearAlgebra.inv(J::HermitianJonesMatrix)
    1/det(J) * HermitianJonesMatrix(J.yy, -J.xy, J.xx)
end

function LinearAlgebra.kron(J1::JonesMatrix, J2::JonesMatrix)
    @SMatrix [J1.xx*J2.xx J1.xx*J2.xy J1.xy*J2.xx J1.xy*J2.xy;
              J1.xx*J2.yx J1.xx*J2.yy J1.xy*J2.yx J1.xy*J2.yy;
              J1.yx*J2.xx J1.yx*J2.xy J1.yy*J2.xx J1.yy*J2.xy;
              J1.yx*J2.yx J1.yx*J2.yy J1.yy*J2.yx J1.yy*J2.yy]
end

@doc doc"""
    congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix)

Compute the congruence transformation of $K$ with respect to $J$. Using this function instead of
explicitly computing $J*K*J'$ guarantees that the final result is exactly Hermitian.

```math
    K \rightarrow JKJ^*
```
"""
function congruence_transform(J::JonesMatrix, K::HermitianJonesMatrix)
    HermitianJonesMatrix(abs2(J.xx)*K.xx + 2real(J.xx*J.xy'*K.xy) + abs2(J.xy)*K.yy,
                         J.xx*J.yx'*K.xx + J.xx*J.yy'*K.xy + J.xy*J.yx'*K.xy' + J.xy*J.yy'*K.yy,
                         abs2(J.yx)*K.xx + 2real(J.yx*J.yy'*K.xy) + abs2(J.yy)*K.yy)
end

function make_hermitian(J::JonesMatrix)
    HermitianJonesMatrix(real(J.xx), 0.5*(J.xy+conj(J.yx)), real(J.yy))
end

function make_hermitian(J::DiagonalJonesMatrix)
    HermitianJonesMatrix(real(J.xx), 0, real(J.yy))
end

make_hermitian(J::HermitianJonesMatrix) = J

