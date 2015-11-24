type Observation
    # Observation lat/y
    y::Float64
    # Observation lon/x
    x::Float64
    # Observation z
    z::Float64
    # Observation time
    t::DateTime
    # Observation kind
    obkind::ASCIIString

    # Value
    value::Float64
    # Error variance
    error::Float64

    # What to do with this ob
    # Assimilate this?
    assim_this::Bool
    assimilated::Bool

    # After assimilation, store info about ob
    prior_mean::Float64
    prior_var::Float64
    post_mean::Float64
    post_var::Float64

    # Constructor function
    function Observation(value, error, obkind="STATE_VARIABLE", t=DateTime(1970), z=0.0, y=0.0, x=0.0)
        new(value,error,obkind,t,z,y,x)
    end
end


