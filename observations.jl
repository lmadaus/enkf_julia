type Observation
    # Observation time
    t::DateTime
    # Observation z
    z::Float64
    # Observation lat/y
    y::Float64
    # Observation lon/x
    x::Float64


    # Observation kind
    obkind::AbstractString

    # Value
    value::Float64
    # Error variance
    error::Float64

    # What to do with this ob
    # Assimilate this?
    assim_this::Bool
    assimilated::Bool

    # After assimilation, store info about ob
    prior_mean::Nullable{Float64}
    prior_var::Nullable{Float64}
    post_mean::Nullable{Float64}
    post_var::Nullable{Float64}

    # Constructor function
    function Observation(value, error, obkind, time, z, y, x)
        # Correct longitude
        if x < 0.0
          x = 360.0 + x
        end
        new(time,z,y,x,obkind,value,error,false,false)
    end

    # Make this prettier to display
    function Base.show(io::IO, ob::Observation)
      print(io, "$(ob.t) $(ob.y) $(ob.x) $(ob.z) $(ob.obkind) $(ob.value)")
    end

end


