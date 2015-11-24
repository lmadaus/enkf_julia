include("kind_definitions.jl")
using NearestNeighbors


function interp_4d(meta_tree, state, ob)
  #=
  The full 4d interpolation of the state to the ob will go here

  Takes the meta_tree, state array and an observation type as
  input.  Uses the meta_tree to find the 10 nearest neighbor points
  (which should be in t,z,y,x; assumes kind is in the ob).
  Should be 4 in x-y, 2 in z, 2 in t.  Interpolates based on
  Euclidean distance

  RETURNS: a (1xNens) array of the ensemble estimate of the ob
  =#
  # Define the observation's coordinates

  # For time, subtract from epoch
  if typeof(ob.t) == DateTime
    obtime = Int(ob.t - DateTime(1970)) / 1000 # From milliseconds to seconds
  end
  println(obtime)
  ob_point = [kind_indices[ob.obkind], obtime, ob.z, ob.y, ob.x]

  # Find indices and distances to 10 nearest neighbors in the tree
  idxs, dists = knn(meta_tree, ob_point, 10, true)
  println(idxs)
  println(dists)
  # Find sum of distances for weights
  tot_dist = sum(dists)
  weights = dists / tot_dist

  # Now get the ensemble state at the indexes given by idxs and
  # multiply by the weights
  state_estimate = zeros(size(state)[2])
  for (idx,w) in zip(idxs,weights)
    weighted = state[idx,:] * w
    state_estimate += weighted'
  end

  return state_estimate

end

function build_state_meta_tree(state_meta)
  #=
  Function that uses the NearestNeighbors package to build
  a tree for efficiently locating the surrounding points of
  a given type,t,z,y,x
  =#
  kdtree = KDTree(state_meta'; leafsize=10)
end


function interp_time(target::DateTime, val1::Float64, time1::DateTime, val2::Float64, time2::DateTime)
  #=
  Given values at two different times (val1, time1) and (val2, time2), do a linear
  interpolation to Datetime 'target'
  =#

  # First be sure that target is between time1 and time2
  if (target < time1) || (target > time2)
    println("ERROR interp_time: target date $target is not between $time1 and $time2")
    return
  end

  # Calculate our linear fit
  dVal = val2 - val1
  dTime = Int(time2 - time1)
  dTarget = Int(target - time1)
  out = val1 + (dVal/dTime) * dTarget

end



function haversine(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)
  #=
  Uses the Haversine formula to calculate the distance between two
  lat/lon pairs (in degrees) in kilometers assuming spherical earth
  =#
  out :: Float64
  R = 6372.8 # Radius of the Earth in km
  out = 2 * R * asin(sqrt(sind((lat2-lat1)/2)^2 + cosd(lat1) * cosd(lat2) * sind((lon2-lon1)/2)^2))

end


function localize(ob::Observation, state_meta, cutoff=100.0, method="GC")
  #=
  Return an array the length of the state with corresponding
  weights to each variable based on distance from the observation

  Requires:
    ob -> Observation object
    state_meta -> State metadata array (Nstate x 5)
    cutoff -> Localization cutoff distance in km (this is distance
              where the weights go to zero)
    method -> Localization method.  Current options are:
              "GC" -- Gaspari-Cohn

  =#

  # Build an array for our distances
  Nstate = size(state_meta)[1]
  localization = zeros(Float64, Nstate)

  # Loop through the state meta data and compute the distance
  # from ob to each point (2-d for now)
  for n=1:Nstate
    dist = haversine(ob.y,ob.x,state_meta[n,4],state_meta[n,5])
    r = dist / (cutoff * 0.5)
    if method == "GC"
      if r <= 1.0
        localization[n] = (((-0.25*r+0.5)*r+0.625)*r-5.0/3.0) * r^2 + 1.0
      elseif (r > 1.0) & (r < 2.0)
        localization[n] = ((((r/12.0 - 0.5)*r + 0.625) * r +
                          5.0/3.0)*r-5.0)*r + 4.0 - 2.0 /
                          (3.0 * r)
      end # GC location check


    end #method

  end # end for n=[1:Nstate]

  return localization


end


