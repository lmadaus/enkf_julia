type Observation
    # Observation lat/y
    y::Float64
    # Observation lon/x
    x::Float64
    # Observation z
    z::Float64
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
    function Observation(y,x,value,error,z=0.0)
        new(y,x,z,value,error)
    end
end



function enkf_update(prior_state, obs; inflate=1.0, loc=false)
    #= 
    Master EnKF update algorithm using EnSRF form of equations

    Requires:
        prior_state -> An array of Nstate x Nmems that is the
                       ensemble state to update. This code
                       currently ASSUMES the prior estimates 
                       of the observations are pre-computed and
                       all appended to the end of the ensemble
                       in the same order they are listed in 
                       the obs variable
        obs         -> A list of observation objects
                       (See the Observation type definition
                       in this file for how to specify them)
    Optional:
        inflate     -> For now, a simple numerical value
                       that will be multiplied by the covariance
                       matrix in the Kalman gain to inflate the
                       spread.  Default is 1.0 (no inflation)
        loc         -> Boolean true/false that will eventually
                       be changed to something more meaningful
                       once localization is figured out.  For
                       now it is not possible to localize.
        
    Returns:
        A tuple of (post_state, obs) with the posterior state
        and the observations used


    =#

    # Figure out the ensemble size
    Nens = size(prior_state)[2]
    Nobs = size(obs)[1]
    # We assume ob estimates are appended to state
    Nstate = size(prior_state)[1] - Nobs
    println("Ensemble size: ", Nens, " State vector size: ", Nstate)

    # Get the ensemble mean and perturbations
    # xbm is background mean, Xbp is background perts
    xbm = mean(prior_state,2)
    Xbp = broadcast(-, prior_state, xbm)
    # Copy these to the post state
    # xam is analysis mean, Xap is analysis perts
    xam = deepcopy(xbm)
    Xap = deepcopy(Xbp)

    # Loop through observations
    for obvs in enumerate(obs)
       ob = obvs[2]
       obnum=obvs[1]
       # Reset the post to be prior
       xbm = xam
       Xbp = Xap
   
       # We construct our H matrix knowing which ob
       # number we are on
       H = zeros(Nstate+Nobs)
       H[Nstate+obnum] = 1
       H = H'
       # Find Ye (HXb, prior estimate of ob)
       yem = H * xbm
       yep = H * Xbp


       # Compute innovation
       innov = ob.value - yem
     
       # check to see if we assimilate this
       # If not, skip the ob
       if !ob.assim_this
           continue
       end


       # Now the innovation variance (denominator of gain)
       kdenom = ob.error + var(yep)

       # Numerator of gain is covariance between
       # ensemble members and estimate of ob
       kcov = (Xbp * yep') / (Nens - 1.0)

       # Option to inflate the covariances by a certain factor
       if inflate != 1.0
           kcov *= inflate
       end

       # Localization here (to be added)

       # Compute Kalman gain
       K = kcov / kdenom

       # Update the mean
       xam = xbm + K * innov

       # Update member perturbations
       # With square root factor
       beta = 1.0 / (1.0 + sqrt(ob.error/(var(yep)+ob.error)))
       K *= beta

       # Update the perturbations
       Xap = Xbp - K * yep

       # record the ob as assimilated
       ob.assimilated = true

       
    end

    # Return the updated state and obs
    post_state = broadcast(+,xam,Xap)

    return (post_state, obs)


end
