using H3.API
using CSV
using Printf
using TOML
using DataFrames


# Define compatibility shim if `geoToH3` doesn't exist, but `latLngToCell` does
if !(@isdefined geoToH3) && (@isdefined latLngToCell)
    @info "Defining shim: geoToH3 from latLngToCell"
    function geoToH3(lat::Real, lng::Real, res::Integer)
        return latLngToCell((lat, lng), res)
    end
end


function get_gr_params(config::Dict, source_id::String)
"""
    get_gr_params(config, source_id

"""
    if haskey(config["sources"], source_id)
        if haskey(config["sources"][source_id], "agr")
            agr = config["sources"][source_id]["agr"]
        else
            msg = @sprintf("Source %s does not have the agr attribute", source_id)
            throw(error(msg))
        end
        if haskey(config["sources"][source_id], "bgr")
            bgr = config["sources"][source_id]["bgr"]
        else
            msg = @sprintf("Source %s does not have the bgr attribute", source_id)
            throw(error(msg))
        end
    end
    return agr, bgr
end


function get_rate_double_truncated_gr(agr::Number, bgr::Number, mag::Number, binw::Number=0.1)
"""
    get_rate_double_truncated_gr(agr, bgr, mag, binw)
    
Using a double truncated GR relationship, computes the rate of occurrence of 
a given magnitude.

# Examples
```julia-repl
julia> get_rate_double_truncated_gr(4.0, 1.0, 5.0, 0.1)
```
"""
    lowm = floor(mag/binw)*binw
    uppm = lowm + binw
    lorate = 10^(agr-bgr*lowm)
    uprate = 10^(agr-bgr*uppm)
    return lorate-uprate
end


function get_completeness_table(config::Dict, source_id::String)
"""
    get_completeness_table(config, fname_config[, fname_out])

Given a dictionary `config` - containing configuration information - and a 
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
                     max_year::Number=3000, weighting::String="one")
"""
boxcounting(fname, hres, fname_h3_to_zone, fname_coonfig, folder_out, max_year, weighting)

Using H3 as a reference, counts the number of earthquakes in each cell.
"""

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

        # This is the default weight assigned to an event
        if weighting == "one"
            
            rate = 1.0

        elseif weighting == "completeness"

            # Check the existance of the configuration file
            @assert length(fname_config) > 0 

            # Check if the mapping between h3 and zones contains the ID for 
            # the current cell and get the corresponding completeness table
            if haskey(mapping, base)
                key = mapping[base]
                compl = get_completeness_table(config, key) 
            else
                continue
                # compl = get_completeness_table(config, "-0")
            end

            # Get the rate
            rate = get_rate_from_completeness(compl, coo[4], coo[3], max_year)
            
        elseif weighting == "mfd"

            # Check the existance of the configuration file
            @assert length(fname_config) > 0 

            # Check if the mapping between h3 and zones contains the ID for 
            # the current cell and get the corresponding completeness table
            if haskey(mapping, base)
                key = mapping[base]
                agr, bgr = get_gr_params(config, key) 
            else
                msg = @sprintf("The h3 to src mapping does not contain this ID: %d", base)
                println(msg)
                continue
            end

            # Get the rate
            binw = 0.1
            rate = get_rate_double_truncated_gr(agr, bgr, coo[4], binw)
        else
            msg = @sprintf("Unknown weighting option: %s", weighting) 
            throw(error(msg))
        end

        # Check the rate assigned to the current earthquake
        if rate < 1e-20
            continue
        end
        
        # Update the counting dictionary
        if haskey(count, base)
            count[base] += rate
        else
            count[base] = rate
        end
    end

    # Creating the putput dataframes
    coo = DataFrame(h3idx=UInt64[], lon=Float64[], lat=Float64[], count=Float64[])
    for k in keys(count)
        geo = h3ToGeo(k)
        push!(coo, [k, rad2deg(geo.lon), rad2deg(geo.lat), count[k]])
    end

    # Write output .csv files
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

    @assert abs(1.0 - sum(weights)) < 1e-10
    
    # Total activity rate
    total_activity_rate = 10^aGR
    
    # Computing local aGR
    aGR_points = log10.(total_activity_rate*weights)
    bGR_points = ones(size(aGR_points))*bGR
    
    # Creating output file
    outdf = DataFrame(lon=points_df.lon, lat=points_df.lat, agr=aGR_points, 
                      bgr=bGR_points)
    CSV.write(fname_out, outdf);
    
end
