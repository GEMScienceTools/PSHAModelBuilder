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
        idxs = kRing(base, maxdistk)

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
                println("wrong index ",base," ",maxdistk)
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

            if lats[idx[2][1]] > 60
                println(idx[2][1], " ", lons[idx[2][1]], " ",lats[idx[2][1]], " ",maxdistk)
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


function boxcounting(fname::String, h3res::Int)
"""
    
    boxcounting(fname, h3res) 
    
Using the H3 library it counts the earthquakes from the catalogue included
in the .csv file defined by `fname`. 

#Examples
```julia-repl
julia> boxcounting('catalogue.csv', 5)
```
"""

    println("Edge length ", edgeLengthKm(h3res))

    fname_out = joinpath("./tmp", @sprintf("box_counting_%s", basename(fname)))
    fname_out_h3 = joinpath("./tmp", @sprintf("box_counting_h3_%s", basename(fname)))

    # Load the catalogue 
    df = DataFrame(CSV.File(fname));

    count = Dict{UInt64,Int64}()
    for coo in zip(df.longitude, df.latitude)

        # GeoCoord takes a lat and a lon
        base = geoToH3(GeoCoord(deg2rad(coo[2]), deg2rad(coo[1])), h3res);
        
        # Update the counting
        if haskey(count, base)
            count[base] += 1
        else
            count[base] = 1
        end
    end

    coo = DataFrame(h3idx=UInt64[], lon=Float64[], lat=Float64[], count = Int32[])
    for k in keys(count)
        geo = h3ToGeo(k)
        push!(coo, [k, rad2deg(geo.lon), rad2deg(geo.lat), count[k]])
    end

    # Write output .csv files
    CSV.write(fname_out, select(coo, :lon, :lat, :count));
    CSV.write(fname_out_h3, select(coo, :h3idx, :lon, :lat, :count));

    println("Writing ", fname_out)
    println("Writing ", fname_out_h3)

end


function distribute_total_rates(aGR::Float64, bGR::Float64, 
                                fname_points::String, fname_out::String)
"""
    
    distribute_total_rates(aGR, bGR, fname_points, fname_out)

Distributes the seismicity specified by the `aGR` and `bGR` parameters
over an irregular grid of `fname_points`. `fname_points` is a .csv file 
with three columns: lon, lat, nocc
The output goes into the file `fname_out`, a .csv formatted file with 
the following columns: lon, lat, aGR, bGR

# Examples
```julia-repl
julia> distribute_total_rates(4.0, 1.0, './smooth.csv'., 'sources.xml')
```
"""
    
    # Load the points
    points_df = DataFrame(CSV.File(fname_points));
    
    # Normalize the weights;
    normalizing_factor = sum(points_df.nocc)
    weights = points_df.nocc ./ normalizing_factor

    # Checking that weights are okay
    @assert abs(sum(weights)-1.0) < 1e-7 "Weights does not sum to unity"
    
    # Total activity rate
    total_activity_rate = 10^aGR
    
    # Computing local aGR
    aGR_points = log10.(total_activity_rate*weights)
    bGR_points = ones(size(aGR_points)) * bGR
    
    # Creating output file
    outdf = DataFrame(lon=points_df.lon, lat=points_df.lat, agr=aGR_points, 
                      bgr=bGR_points)
    CSV.write(fname_out, outdf);
    
end
