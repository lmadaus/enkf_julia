using NetCDF
include("kind_definitions.jl")

kind_names = Dict{ASCIIString, ASCIIString}(
  "TEMPERATURE_2M" => "TMP_2maboveground",


  )


function build_state(vars, model_path=".")
  #=
  Given a list of variables, grab these from all ensemble members
  Build the state meta data to accompany it
  =#

  # Get all ensemble members
  files = searchdir(model_path, ".nc")
  Nens = length(files)
  state_meta =[]
  cur_state = []
  full_state = []
  varsize = 0
  # Loop though each file and variable
  for (fnum,f) in enumerate(files)
    # Loop through each variable
    for (varnum, var) in enumerate(vars)
      var_index = kind_indices[var]
      values = ncread(f, kind_names[var])

      # If this is the first ensemble member, prepare the metadata
      if fnum == 1
        # Build the dimensions for this variable
        # Start with variable index
        varsize = size(values)
        varkindarr = fill(var_index, varsize)
        varkindarr = reshape(varkindarr, prod(varsize))
        # Now we need times
        times = ncread(f,"time")
        timearr = broadcast!(identity, ones(varsize), times[1,1,:])
        #println(size(timearr))
        #println(varsize)

        timearr = reshape(timearr, prod(varsize))
        # For SREF z is zero everywhere for now
        zarr = zeros(varsize)
        zarr = reshape(zarr, prod(varsize))
        # Latitudes
        lats = ncread(f,"latitude")
        latarr = broadcast(*, lats, ones(varsize))
        latarr = reshape(latarr, prod(varsize))
        # Longitudes
        lons = ncread(f,"longitude")
        lonarr = broadcast(*, lons, ones(varsize))
        lonarr = reshape(lonarr, prod(varsize))



        # Now we need to concatenate these
        this_block = hcat(varkindarr, timearr, zarr, latarr, lonarr)

        # If this is the first variable, this is now the metadata
        # Otherwise, append to the existing metadata
        if varnum == 1
          state_meta = this_block
        else
          state_meta = vcat(state_meta, this_block)
        end # if varnum == 1
      end # End first ensemble member

      # Build the state for this member
      # Reshape actual values here
      values = reshape(values, prod(varsize))
      if varnum == 1
          cur_state = values
      else
          cur_state = vcat(cur_state, values)
      end # if varnum == 1
      #println("Variable ", var," ", varnum)
      #println(size(cur_state))



    end # end variable loop
    println("State meta shape ", size(state_meta))
    # Now if this is the first ensemble member, Broaden out the entire state array
    if fnum == 1
      Nstate = length(cur_state)
      full_state = zeros((Nstate, Nens))
    end # if fnum == 1

    full_state[:,fnum] = cur_state

  end # end file loop

  # Return the full_state and state_meta
  return (full_state, state_meta)


end

searchdir(path,key) = filter(x->contains(x,key), readdir(path))
