using Plots
using DelimitedFiles
using Measures

const USAGE = """
USAGE: julia plotwfs.jl [OPTIONS] potential.dat energies.dat wfs.dat

  Plots the potential energy curve V(R) and each wavefunction shifted up by its energy (real part).

OPTIONS
  --runits U      length units of the potential file: bohr, au, angstrom               [bohr]
  --eunits U      energy units of the potential file: hartree, au, ev, invcm, rydberg  [hartree]
  --states STATES which states to plot, by vibrational quantum number v                [all]
  --scale  scl    wavefunction scale s, in energy units per unit of ψ                  [auto]
  --density       plot |ψ|² instead of ψ
  --vmax E        upper limit of the energy axis                                       [highest E plotted + a bit]
  --rmin R        lower limit of the R axis, in the units of the wavefunction file     [min in file]
  --rmax R        upper limit of the R axis, in the units of the wavefunction file     [max in file]
  -o file         output file                                                          [wfs.png]
  -h, --help      show this helpful message (:

  With complex wavefunctions, by default the real parts of ψ and E are plotted
"""

const LENGTH_TO_BOHR = Dict(
    "bohr"     => 1.0
  , "au"       => 1.0
  , "angstrom" => 1 / 0.529177210903
 )
const ENERGY_TO_HARTREE = Dict(
    "hartree" => 1.0
  , "au"      => 1.0
  , "ev"      => 1 / 27.211386245988
  , "invcm"   => 1 / 219474.6313632
  , "rydberg" => 0.5
)
const UNIT_LABELS = Dict(
    "ev"       => "eV"
  , "invcm"    => "cm⁻¹"
  , "hartree"  => "hartree"
  , "au"       => "au"
  , "rydberg"  => "Ry"
  , "bohr"     => "bohr"
  , "angstrom" => "Å"
 )

"Returns the first match of `regex` in the `#` header lines of `path` lowercased. Returns `nothing` if no match"
function header_value(path::AbstractString, regex::Regex)
  for line in eachline(path)
    startswith(line, "#") || break
    m = match(regex, line)
    m === nothing || return lowercase(m.captures[1])
  end
  return nothing
end

function unit_factor(table, unit, who)
  haskey(table, unit) || error("Unknown $who unit '$unit'. Try one of: $(join(sort(collect(keys(table))), ", "))")
  return table[unit]
end

readnumbers(path) = Float64.(readdlm(path; comments=true, comment_char='#'))

"Energies file -> (real energies, units, is_complex)"
function read_energies(path)
  data  = readnumbers(path)
  ncols = size(data, 2)
  ncols in (3, 5) || error("$path does not look like a SEECS energies file (expected 3 or 5 columns, got $ncols)")
  units = something(header_value(path, r"energy units:\s*(\S+)"), "hartree")
  return data[:, 2], units, ncols == 5
end

"Wavefunctions file -> (R [nR], ψ [nR × nstates], length units)"
function read_wfs(path, is_complex)
  data  = readnumbers(path)
  R     = data[:, 1]
  cols  = data[:, 2:end]
  if is_complex
    iseven(size(cols, 2)) || error("$path: expected Re/Im column pairs for complex wavefunctions, but got an odd number ($(size(cols, 2))) of columns")
    psi = cols[:, 1:2:end] .+ im .* cols[:, 2:2:end]
  else
    psi = cols
  end
  units = something(header_value(path, r"distance units:\s*(\S+)"), "bohr")
  return R, psi, units
end

function parse_states(s)
  bad() = error("--states must look like 3, 0:10, 0:5:20 or 0,3,7, got '$s'")
  s = replace(s, " " => "")
  occursin(",", s) && return sort(unique(reduce(vcat, parse_states.(split(s, ","; keepempty=false)))))
  parts = split(s, ":")
  all(p -> occursin(r"^\d+$", p), parts) || bad()
  n = parse.(Int, parts)
  length(n) == 1 && return n
  length(n) == 2 && return collect(n[1]:n[2])
  length(n) == 3 && n[2] > 0 && return collect(n[1]:n[2]:n[3])
  bad()
end

"""
    plot_wavefunctions(Vfile, efile, wfsfile
      ; pot_runits = "bohr"
      , pot_eunits = "hartree"
      , states     = nothing
      , scale      = nothing
      , density    = false
      , vmax       = nothing
      , rmin       = nothing
      , rmax       = nothing
      , outfile    = "wfs.png"
      , kwargs...
      )

Plot the wavefunction(s) and potential at their energies. kwargs are passed to `plot`.
"""
function plot_wavefunctions(
      potfile
    , efile
    , wfsfile
    ; pot_runits = "bohr"
    , pot_eunits = "hartree"
    , states     = nothing
    , scale      = nothing
    , density    = false
    , vmax       = nothing
    , rmin       = nothing
    , rmax       = nothing
    , outfile    = nothing
    , palette    = nothing
    , kwargs...
  )
  E, eunits, is_complex = read_energies(efile)
  R, psi, runits        = read_wfs(wfsfile, is_complex)
  size(psi, 2) == length(E) ||
    error("$wfsfile has $(size(psi, 2)) states but $efile has $(length(E)). Are they from the same run ?")

  # -- the potential, converted to the units of the SEECS output files
  pot = readnumbers(potfile)
  size(pot, 2) >= 2 || error("$potfile must have at least two columns: R V(R)")
  Rpot = pot[:, 1] .* unit_factor(LENGTH_TO_BOHR, lowercase(pot_runits), "length") ./ unit_factor(LENGTH_TO_BOHR, runits, "length")
  Vpot = pot[:, 2] .* unit_factor(ENERGY_TO_HARTREE, lowercase(pot_eunits), "energy") ./ unit_factor(ENERGY_TO_HARTREE, eunits, "energy")

  # -- which states; v is 0-based, columns are 1-based
  vs = states === nothing ? collect(0:length(E)-1) : collect(states)
  isempty(vs) && error("No states selected")
  (minimum(vs) >= 0 && maximum(vs) <= length(E) - 1) || error("States $vs out of range: the files hold v = 0 to $(length(E)-1)")
  cols = vs .+ 1

  shape = density ? abs2.(psi[:, cols]) : real.(psi[:, cols])

  # -- automatic scale: the largest plotted wavefunction spans 0.4 of the mean spacing of the plotted levels
  if scale === nothing
    spacing = length(cols) > 1 ? (E[cols[end]] - E[cols[1]]) / (length(cols) - 1) : 0.1 * abs(E[cols[1]] - minimum(Vpot))
    peak    = maximum(abs, shape)
    scale   = peak > 0 ? 0.4 * spacing / peak : 1.0
  end

  # -- R range: where the plotted states are not negligible (above 1e-3 of the largest), plus a margin
  significant = vec(any(abs.(shape) .> 1e-3 * maximum(abs, shape); dims=2))
  Rlo, Rhi    = extrema(R[significant])
  margin      = 0.1 * (Rhi - Rlo)
  rmin = rmin === nothing ? max(minimum(R), Rlo - margin) : rmin
  rmax = rmax === nothing ? min(maximum(R), Rhi + margin) : rmax

  # -- E range: from the well bottom to above the highest plotted level
  inrange = (Rpot .>= rmin) .& (Rpot .<= rmax)
  any(inrange) || error("The potential file has no points between R = $rmin and $rmax $runits")
  vmin    = minimum(Vpot[inrange])
  top     = maximum(E[cols]) + maximum(abs, scale .* shape)
  vmax    = vmax === nothing ? top + 0.05 * (top - vmin) : vmax

  rl = get(UNIT_LABELS, runits, runits)
  el = get(UNIT_LABELS, eunits, eunits)
  p = plot(
      Rpot, Vpot
    ; label      = "V(R)"
    , color      = :black
    , linewidth  = 2
    , xlims      = (rmin, rmax)
    , ylims      = (vmin - 0.02 * (vmax - vmin), vmax)
    , xlabel     = "R [$rl]"
    , ylabel     = "E [$el]"
    , title      = density ? "|ψ_v(R)|²" : (is_complex ? "Re ψ_v(R) (ECS)" : "ψ_v(R)")
    , legend     = :best
    , framestyle = :box
    , left_margin   = 8mm
    , bottom_margin = 6mm
    , right_margin  = 4mm
    , top_margin    = 4mm
    , size       = (1100, 750)
    , dpi        = 300
  )
  palette = palette !== nothing ? palette : cgrad(:viridis, max(length(cols), 2); categorical=true)
  for (n, c) in enumerate(cols)
    hline!(p, [E[c]]; label=false, color=palette[n], linestyle=:dot, alpha=0.6)
    plot!(p, R, E[c] .+ scale .* shape[:, n]; label="v = $(c-1)", color=palette[n], kwargs...)
  end
  if outfile !== nothing
    savefig(p, outfile)
    println("Wrote plot to $outfile (wavefunction scale = $(round(scale; sigdigits=3)) $el per unit of ψ$(density ? "²" : ""))")
  end
  return p
end
