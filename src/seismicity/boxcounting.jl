using H3.API
using CSV
using Printf
using TOML
using DataFrames

function get_completeness_table(config::Dict, source_id::String)
"""
    get_completeness_table(config, fname_config[, fname_out])

Given a dictionary `config` containing configuration information and a 
source ID, this function provides the corresponding completeness table. 
When a source does not have a completeness table, we use the default one.

# Examples
```julia-repl
julia> get_completeness_table(config, "1")
```
"""
    if (haskey(config["sources"], source_id) &&
		haskey(config["sources"][source_id], "completeness_table"))
        compl = config["sources"][source_id]["completeness_table"]
    else
        compl = config["default"]["completeness_table"]
    end
    return compl
end

function get_rate_from_completeness(compl::Array, mag::Number, year::Number, year_end::Number)
    i = size(compl)[1]
    while i > 0 
        if mag >= compl[i][2]
            if year >= compl[i][1]
                return 1.0/(year_end-compl[i][1])
            end
        end
        i -= 1
    end
    return 0.0
end


function boxcounting(fname::String, h3res::Int, fname_h3_to_zone::String="",
                     fname_config::String="", folder_out::String="", 
					 max_year::Number=3000)

    fname_out = joinpath(folder_out, @sprintf("box_counting_%s", basename(fname)))
    fname_out_h3 = joinpath(folder_out, @sprintf("box_counting_h3_%s", basename(fname)))

    # Read the file containing the mapping between h3 cells and zones. This is 
    # used to select information about completeness
    mapping = Dict{UInt64,String}()
    if length(fname_h3_to_zone) > 0
        open(fname_h3_to_zone) do f
            for row in eachline(f)
                tmp = split(row, ",")
                mapping[parse(UInt64, tmp[1], base=16)] = tmp[2]
            end
        end
    end
  
    # Read the TOML configuration file
    if length(fname_config) > 0
        @assert length(fname_h3_to_zone) > 0
        config = TOML.parsefile(fname_config)
    end

    # Load the .csv file containing the catalogue 
    df = DataFrame(CSV.File(fname));

    # Loop over the earthquakes in the catalogue
    count = Dict{UInt64,Float64}()
    for coo in zip(df.longitude, df.latitude, df.year, df.magnitude)

        # GeoCoord takes a lat and a lon
        base = geoToH3(GeoCoord(deg2rad(coo[2]), deg2rad(coo[1])), h3res);

        rate = 1.0
        if length(fname_config) > 0 
			if haskey(mapping, base)
            	compl = get_completeness_table(config, mapping[base]) 
			else
				compl = get_completeness_table(config, "-0")
        	end
			rate = get_rate_from_completeness(compl, coo[4], coo[3], max_year)
		end

        if rate < 1e-20
            continue
        end
        
        # Update the counting
        if haskey(count, base)
            count[base] += rate
        else
            count[base] = rate
        end
    end

    coo = DataFrame(h3idx=UInt64[], lon=Float64[], lat=Float64[], count=Float64[])
    for k in keys(count)
        geo = h3ToGeo(k)
        push!(coo, [k, rad2deg(geo.lon), rad2deg(geo.lat), count[k]])
    end

    # Write output .csv files
    println(folder_out)
    if length(folder_out) > 0
        CSV.write(fname_out, select(coo, :lon, :lat, :count));
        CSV.write(fname_out_h3, select(coo, :h3idx, :lon, :lat, :count));
        println("Writing ", fname_out)
        println("Writing ", fname_out_h3)
    end

    return coo

end


function distribute_total_rates(aGR::Float64, bGR::Float64, fname_in::String, fname_out::String)
    """
    
    distribute_total_rates(aGR, bGR, fname_in, sources)

    Distributes the seismicity specified by the `aGR` and `bGR` parameters
    over an irregular grid of `fname_in`. `fname_in` is a .csv file with three
    columns: lon, lat, nocc
    The output goes into the file `fname_out`, a .csv formatted file with 
    the following columns: lon, lat, aGR, bGR

    # Examples
    ```julia-repl
    julia> distribute_total_rates(4.0, 1.0, './smooth.csv'., 'sources.xml')
    ```
    """
    
    # Load the points
    points_df = DataFrame(CSV.File(fname_in));
    
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
