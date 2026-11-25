using NumericalIntegration: integrate
using DelimitedFiles: writedlm, readdlm

###############################################################################
#                                                                             #
#    state-selected thermal convolution of raw cross sections behaving as 1/E #
#    at threshold with a Maxwell-Boltzmann distribution.                      #
#    See, e.g.,                                                               #
#                                                                             #
#        https://doi.org/10.1063/1.2784275                                    #
#                                                                             #
#    reads a data file formatted as :                                         #
#                                                                             #
#      x  y                                                                   #
#      .  .                                                                   #
#      .  .                                                                   #
#      .  .                                                                   #
#                                                                             #
#    where x is an energy and y is a cross section.                           #
#    Should be used with the wrapper script by the same name.                 #
#                                                                             #
###############################################################################

PROGRAM = basename(@__FILE__)

ln(x) = log(x)

logrange(a,b,n) = (10^x for x in range(log10(a), log10(b), length=n))

au2cm = 0.5291772083e-8 # <---- atomic units -> centimeters (length)
au2s  = 2.4188843e-17 # <------ atomic units -> seconds (time)
α_au2SI = au2cm^3 / au2s # <-- atomic units -> centimeters^3 per second (rate coefficient)
au2K = 3.1577465e5 # <--------- atomic units -> Kelvin (temperature)

function usage(::Int)
  println()
  println(PROGRAM,": Convolves data from a datafile with a function ")
  println()
  println("    usage: ",PROGRAM," [FILE]")
  println()
  return true
end

function main()

  narg = 12

  if(size(ARGS,1) != narg)
    println("STOP: ",PROGRAM," requires ",narg," arguments, but ",string(size(ARGS,1))," were supplied.")
    usage(1)
    exit()
  end

  datafile     = ARGS[1]
  outfile      = ARGS[2]
  Ti           = parse(Float64,ARGS[3])  # <--- lowest  temperature for kinetics data
  Tf           = parse(Float64,ARGS[4])  # <--- highest temperature for kinetics data
  EelMin       = parse(Float64,ARGS[5])  # <--- lowest electron energy to which cross sections are
                                         # | extrapolated for the integral over electron energies.
                                         # | Useful for integrating near 0 where the cross sections
                                         # | have a simple 1/E behavior.
  numExtrap    = parse(Int,ARGS[6])      # <------- number of extrapolated grid points
  nT           = parse(Int,ARGS[7])      # <------- number of grid points
  logx         = parse(Bool,ARGS[8])     # <------ T: logarithmic convolution grid
  extrap       = parse(Bool,ARGS[9])     # <------ T: logarithmic convolution grid
  xsConversion = parse(Float64,ARGS[10]) # <-- : conversion factor supplied by wrapper to get energy in hartree atomic units
  EConversion  = parse(Float64,ARGS[11]) # <-- : conversion factor supplied by wrapper to get cross sections in bohr^2
  comment_char = only(ARGS[12])          # <-- : convert string to character

  # -- check if input datafile exists
  if( ! isfile(datafile) )
    usage(1)
    println("STOP: Need a data file to convole")
    exit()
  end

  # -- assumes x y format, one y column (for now)
  data  = readdlm(datafile, comments = true, comment_char = comment_char)
  xdata = data[:,1] * EConversion
  ydata = data[:,2] * xsConversion ^ 2

  # -- sort electron energies by size if not done so already
  indices = sortperm(xdata)
  EelGrid = xdata[indices]

  # -- get probabilities from cross sections
  probs = ydata[indices] * 2 .* EelGrid / pi

  # -- prepend extrapolated energy grid to approach zero. Assume probabilities are constant at low energies
  Eel0 = first(EelGrid)
  EelGridExtrap = [ nrg for nrg in logrange(EelMin,Eel0,numExtrap)]
  EelGrid = vcat(EelGridExtrap, EelGrid[2:end])
  if( extrap )
    probs = vcat([first(probs) for p in EelGridExtrap[2:end]], probs)
  else
    probs = vcat([0 for p in EelGridExtrap[2:end]], probs)
  end

  # -- define temperature grid
  if(logx)
    kTGrid = [ T for T in logrange(Ti,Tf,nT) ]
  else
    kTGrid = [ T for T in range(Ti,Tf,nT) ]
  end

  numer = [ integrate(EelGrid, probs          .* exp.(- EelGrid ./ kT)) for kT in kTGrid ./ au2K ]
  denom = [ sqrt(pi) / 2 * kT^1.5 for kT in kTGrid ./ au2K ]

  rate = numer ./ denom

  # rate = rate * pi / sqrt(2) * α_au2SI
  rate = rate * pi / sqrt(2) * α_au2SI

  writedlm(outfile, [kTGrid rate])

end

main()
