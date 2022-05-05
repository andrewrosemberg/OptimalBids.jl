```@meta
CurrentModule = OptimalBids
```

```@raw html
<div style="width:100%; height:150px;border-width:4px;border-style:solid;padding-top:25px;
        border-color:#000;border-radius:10px;text-align:center;background-color:#99DDFF;
        color:#000">
    <h3 style="color: black;">Star us on GitHub!</h3>
    <a class="github-button" href="https://github.com/andrewrosemberg/OptimalBids.jl" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star andrewrosemberg/OptimalBids.jl on GitHub" style="margin:auto">Star</a>
    <script async defer src="https://buttons.github.io/buttons.js"></script>
</div>
```

# OptimalBids

Documentation for [OptimalBids](https://github.com/andrewrosemberg/OptimalBids.jl).

Simple Package to help evaluate bids in markets/auctions.

## Installation

```julia
] add OptimalBids
```

## Overview

The framework used is based around the dispatch of the following functions for specified markets (subtypes of the core abstract type `Market`):
`build_market`, `change_bids!`, `clear_market!`, `calculate_profit` (see docstrings for more info - e.g. `@doc change_bids!`). I.e. these are the functions that need to be implemented for each market type.

In addition, the functions `profit_for_bid!`, `profit_curve!` are implemented in a reasonably generic way, allowing users to dispatch them with any new market. (see docstrings for more info - e.g. `@doc profit_for_bid!`).
