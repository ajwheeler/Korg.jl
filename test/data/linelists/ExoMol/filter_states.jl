# Read the transitions file and collect all unique state IDs
function collect_state_ids(trans_file)
    state_ids = Set{Int}()

    open(trans_file) do f
        for line in eachline(f)
            # Skip empty lines
            isempty(strip(line)) && continue

            # Split line on whitespace and get first two columns
            fields = split(line)
            push!(state_ids, parse(Int, fields[1]))
            push!(state_ids, parse(Int, fields[2]))
        end
    end

    return state_ids
end

# Filter states file to only include states referenced in transitions
function filter_states(states_file, output_file, state_ids)
    open(states_file) do input
        open(output_file, "w") do output
            for line in eachline(input)
                # Skip empty lines
                isempty(strip(line)) && continue

                # Try to get state ID from first column
                fields = split(line)
                try
                    state_id = parse(Int, fields[1])
                    # Only write lines where state ID is in our set
                    if state_id in state_ids
                        println(output, line)
                    end
                catch
                    println("Error parsing state ID: $line")
                    # Skip lines that don't start with an integer
                    continue
                end
            end
        end
    end
end

function (@main)(ARGS)
    trans_file = "40Ca-1H__XAB_abridged.trans"
    states_file = "40Ca-1H__XAB.states"
    output_file = "40Ca-1H__XAB_abridged.states"

    state_ids = collect_state_ids(trans_file)
    filter_states(states_file, output_file, state_ids)
end
