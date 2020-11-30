using H3.API
using CSV
using Printf
using DataFrames


function boxcounting(fname::String, h3res::Int)

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


function distribute_total_rates(aGR::Float64, bGR::Float64, points::String, fname_out::String)
    """
    
    distribute_total_rates(aGR, bGR, points, sources)

    Distributes the seismicity specified by the `aGR` and `bGR` parameters
    over an irregular grid of `points`. `points` is a .csv file with three
    columns: lon, lat, nocc
    The output goes into the file `fname_out`, a .csv formatted file with 
    the following columns: lon, lat, aGR, bGR

    # Examples
    ```julia-repl
    julia> distribute_total_rates(4.0, 1.0, './smooth.csv'., 'sources.xml')
    ```
    """
    
    # Load the points
    points_df = DataFrame(CSV.File(fname));
    
    # Normalize the weights;
    normalizing_factor = sum(points_df.nocc)
    weights = points_df.nocc ./ normalizing_factor
    
    # Total activity rate
    total_activity_rate = 10^aGR
    
    # Computing local aGR
    aGR_points = log10.(total_activity_rate*weights)
    bGR_points = ones(size(aGR_points))
    
    # Creating output file
    outdf = DataFrame(lon=points_df.lon, lat=points_df.lat, agr = aGR_points, 
                      bgr = bGR_points)
    CSV.write(fname_out, outdf);
    
end
