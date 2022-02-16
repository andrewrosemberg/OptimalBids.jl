"""
    add_generator(nw_data::Dict, bus_index::Int) -> String

Adds generator at bus with index `bus_index` in grid data (`nw_data`), and returns dictionary key for the new generator.
"""
function add_generator(nw_data::Dict, bus_index::Int)
    # index new generator as the last generator
    gen_idx = length(nw_data["gen"]) + 1

    # error if indexes are not sequential
    !ismissing(get(nw_data["gen"], string(gen_idx), missing)) &&
        error("crazy generator indexes")

    nw_data["gen"][string(gen_idx)] = Dict(
        "ncost" => 1,
        "qc1max" => 0.0,
        "pg" => 0.0,
        "model" => 2,
        "shutdown" => 0.0,
        "startup" => 0.0,
        "qc2max" => 0.0,
        "ramp_agc" => 0.0,
        "qg" => 0.0,
        "gen_bus" => bus_index,
        "pmax" => 0.0,
        "ramp_10" => 0.0,
        "vg" => 0.955,
        "mbase" => 100.0,
        "source_id" => Any["gen", gen_idx],
        "pc2" => 0.0,
        "index" => gen_idx,
        "cost" => 0.0,
        "qmax" => 0.0,
        "gen_status" => 1,
        "qmin" => 0.00,
        "qc1min" => 0.0,
        "qc2min" => 0.0,
        "pc1" => 0.0,
        "ramp_q" => 0.0,
        "ramp_30" => 0.0,
        "pmin" => 0.0,
        "apf" => 0.0,
    )

    return string(gen_idx)
end
