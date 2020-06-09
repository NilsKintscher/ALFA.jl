function Base.show(io::IO, mime::MIME"text/plain", C::Crystal)
    print(io, "Lattice Basis: ")
    show(io, mime, C.L)
    print(io, "\nDomain: ")
    show(io, mime, C.Domain)
    print(io, "\nCodomain: ")
    show(io, mime, C.Codomain)
end

function Base.show(io::IO, mime::MIME"text/plain", m::Multiplier)
    print(io, "Position: ")
    show(io, mime, m.pos)
    print(io, "\nMultiplier: ")
    show(io, mime, m.mat)
end

function Base.show(io::IO, mime::MIME"text/plain", s::SortedSet{Multiplier})
    show(io, mime, collect(s))
end

function Base.show(io::IO, mime::MIME"text/plain", o::CrystalOperator)
    show(io, mime, o.C)
    print(io, "\nMultiplier: ")
    show(io, mime, o.M)
end
