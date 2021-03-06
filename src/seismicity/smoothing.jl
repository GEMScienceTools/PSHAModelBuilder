using H3.API
using CSV
using TOML
using Printf
using DataFrames
using Distributions


function smoothing(fname_count::String, fname_config::String, fname_out::String="")

    model = TOML.parsefile(fname_config)
    maxdistkm = model["kernel_maximum_distance"]
    smoothing_σs = model["kernel_smoothing"]
    smooth = smoothing(fname_count, smoothing_σs, maxdistkm)

    minlo = describe(smooth, :min, cols=:lon).min[1]
    maxla = describe(smooth, :max, cols=:lat).max[1]

    if length(fname_out) == 0
        tmp = split((basename(fname_count)), '.')
        fname_out = joinpath("./tmp", @sprintf("%s_smooth.csv", tmp[1]));  
    end

    CSV.write(fname_out, select(smooth, :lon, :lat, :nocc));
    println("Created         : ", fname_out)
end


function smoothing(fname_count::String, smoothing_σs::Array, maxdistkm::Real)
"""
    smoothing(fname, smoothing_σs, maxdistkm) 
    
The `fname_count` is a .csv file containing for each row a (lon, lat, count)
tuple. This function smooths the count in this file using the gaussian kernels
defined in the array `smoothing_σs`; this array contains N (weight, std) 
tuples where std is a standard deviation [km]. The sum of weights in 
`smoothing_σs` must be equal to 1. The `maxdistkm` is the maximum distance [km]
used to perform the smoothing.

# Examples
```julia-repl
julia> smoothing('count.csv', [[1.0, 20]], 50)
```
"""

    df = DataFrame(CSV.File(fname_count));
    df[!,:h3idx] = convert.(UInt64,df[!,:h3idx]);

    # Find the resolution and according to h3 resolution
    h3res = h3GetResolution(df.h3idx[1])
    println(@sprintf("Edge resolution : %d", h3res))
    println(@sprintf("Edge length     : %.3f km", edgeLengthKm(h3res)))

    maxdistk = Int(ceil(maxdistkm/edgeLengthKm(h3res)))
    println(@sprintf("Max dist k      : %d ", maxdistk))

    nocc = Dict{UInt64,Float32}()
    lons = Dict{UInt64,Float32}()
    lats = Dict{UInt64,Float32}()
    println("Number nodes    : ", length(df.h3idx))

    for tmp in enumerate(zip(df.lon, df.lat, df.count))

        # Cell index to whick the current point (i.e. cell center) belongs to
        base = geoToH3(GeoCoord(deg2rad(tmp[2][2]), deg2rad(tmp[2][1])), h3res);

        # Get the indexes of cells within 'maxdistk' from the cell with index
        # equal to 'base'. 'maxdistk' is the distance in terms of the cell 
        # size
        tmp_idxs = kRing(base, maxdistk)

        # Removing invalid indexes
        idxs = []
        for i in tmp_idxs
            if h3IsValid(i)
                push!(idxs, i)
            end
        end

        dsts = zeros(Float32, length(idxs))
        for idx in enumerate(idxs)
            d = h3Distance(base, idx[2])
            dsts[idx[1]] = d * edgeLengthKm(h3res)
            if dsts[idx[1]] < 1.0
                dsts[idx[1]] = 1.0
            end
        end

        # Gaussian weights
        wei = zeros(size(dsts))
        for smo in smoothing_σs
            wei += pdf.(Normal(0.0, smo[2]), dsts) .* smo[1]
        end

        # Normalising
        wei /= sum(wei)
        sum(wei)-1.0 < 1e-5 || error("weights are wrong") 

        # Updating the array with the smoothing
        for idx in enumerate(zip(idxs, wei))

            if !h3IsValid(idx[2][1])
                println(@sprintf("wrong index obtained from kRing with 0x%08x", base))
            end

            # Updating the nocc, lons and lats dictionaries
            if haskey(nocc, idx[2][1])
                nocc[idx[2][1]] += tmp[2][3] * idx[2][2]
            else
                nocc[idx[2][1]] = tmp[2][3] * idx[2][2]
                geo1 = h3ToGeo(idx[2][1])
                lons[idx[2][1]] = rad2deg(geo1.lon)
                lats[idx[2][1]] = rad2deg(geo1.lat)
            end

        end

    end

    # Save results into file
    smooth = DataFrame(lon=Float64[], lat=Float64[], nocc=Float64[])
    for k in keys(nocc)
        push!(smooth, [lons[k], lats[k], nocc[k]])
    end  
  
    return smooth
end
