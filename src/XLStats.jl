module XLStats

using OrderedCollections
using Statistics
using Parameters
using Printf

export Link, getscore, read_all, filter, point_biserial, indomain

export name, resnames, indexes, ismatch, getlink,
       consistency, deuc, dtop, dmax, nscans, 
       maxscore1, avgscore1,
       maxscore2, avgscore2,
       hasxic, maxxic, avgxic, 
       nspecies, count_nspecies

#
# Convert one letter and three-letter amino acid codes
#
const oneletter = Dict("CYS"=>"C", "ASP"=>"D", "SER"=>"S", "GLN"=>"Q", "LYS"=>"K",
                       "ILE"=>"I", "PRO"=>"P", "THR"=>"T", "PHE"=>"F", "ASN"=>"N", 
                       "GLY"=>"G", "HIS"=>"H", "LEU"=>"L", "ARG"=>"R", "TRP"=>"W", 
                       "ALA"=>"A", "VAL"=>"V", "GLU"=>"E", "TYR"=>"Y", "MET"=>"M")
const threeletter = Dict( v => k for (k,v) in oneletter )

#
# List of possible data of the scans and the associated symbols in the Link struct
# needs to be ordered because it will define the substitution sequence when parsing
# the Scan line, and score1 and score2 conflict.
#
const data_list = OrderedDict( "Scan:"                => :index, 
                               "Peptide:"             => :peptide,
                               "Secondary Score:"     => :score2,
                               "Score:"               => :score1,
                               "XL Position 1:"       => :pos1,   
                               "XL Position 2:"       => :pos2,
                               "Source File:"         => :source_file,
                               "Experimental M+H:"    => :mplush,
                               "Precursor Charge:"    => :prec_charge,
                               "Peaks Matched alpha:" => :matched_alpha,
                               "Peaks Matched beta:"  => :matched_beta,
                               "Retention Time:"      => :retention_time,
                               "Assessment:"          => :assessment )

#
# Structure to contain scan data
#
@with_kw mutable struct Scan
  index :: Int = 0
  pep1 :: String = ""
  pep2 :: String = ""
  score1 :: Float64 = 0.
  score2 :: Float64 = 0.
  pos1 :: Int = 0
  pos2 :: Int = 0
  source_file :: String = ""
  mplush :: Float64 = 0
  prec_charge :: Int = 0
  matched_alpha :: Int = 0
  matched_beta :: Int = 0 
  retention_time :: Float64 = 0
  xic :: Float64 = -1
end

#
# Structure that contain each link data
#
@with_kw mutable struct Link

  name::AbstractString = "None"

  consistency::Bool = false
  deuc::Float64 = 0.
  dtop::Float64 = 0.
  dmax::Float64 = 0.

  nspecies :: Int = 0
  hasxic::Bool = false
  scans :: Vector{Scan} 

end

function Link(name::String,nscans::Int)  
  str = split(name,'-')
  i = parse(Int,str[1][2:end])
  j = parse(Int,str[2][2:end])
  scans = [ Scan() for i in 1:nscans ]
  Link(name=name, scans=scans)
end

# 
# syntax sugars for getfield on array of links
#
name(link::Link) = link.name 
name(links::Vector{Link}) = name.(links)

function resnames(linkname::AbstractString;type=3)
  name = split(linkname,"-")
  if type == 3
    return threeletter[name[1][1:1]], threeletter[name[2][1:1]]
  elseif type == 1
    return name[1][1:1], name[2][1:1]
  else
    error(" type must be 1 or 3 in resnames. ")
  end
end
resnames(link::Link;type=3) = resnames(link.name,type=type)
  
function indexes(linkname::AbstractString)
  name = split(linkname,"-")
  return parse(Int,name[1][2:end]), parse(Int,name[2][2:end])
end
indexes(link::Link) = indexes(link.name)

function ismatch(name1::AbstractString,name2::AbstractString)
   ( ismatch(resnames(name1),resnames(name2)) &&
     ismatch(indexes(name1),indexes(name2)) )
end

function ismatch(x::Tuple{T,T},y::Tuple{T,T}) where T
  ( (x[1] == y[1] && x[2] == y[2]) ||
    (x[2] == y[1] && x[1] == y[2]) )
end

# 
# Allow retrieving links by name (much slower, of course)
#
import Base.getindex
getindex(links::Vector{Link},name::AbstractString) = links[findfirst( link -> ismatch(link.name,name), links )]

consistency(link::Link;tol=0.) = (0. <= link.dtop <= (link.dmax + tol))
consistency(links::Vector{Link};tol=0) = consistency.(links,tol=tol)

deuc(link::Link) = link.deuc
deuc(links::Vector{Link}) = deuc.(links)

dtop(link::Link) = link.dtop
dtop(links::Vector{Link}) = dtop.(links)

dmax(link::Link) = link.dmax
dmax(links::Vector{Link}) = dmax.(links)

nscans(link::Link) = length(link.scans)
nscans(links::Vector{Link}) = nscans.(links)

#
# Get maximum value of a field of the scans of a link
#
function getmax(link::Link,field::Symbol)
  nscans(link) == 0 && return missing
  valmax = getfield(link.scans[1],field)
  for iscan in 2:nscans(link)
     valmax = max(valmax,getfield(link.scans[iscan],field))
  end
  return valmax
end

#
# Get the mean value of a field of the scans of a link
#
function getmean(link::Link,field::Symbol)
  nscans(link) == 0 && return missing
  meanval = zero(typeof(getfield(link.scans[1],field)))
  for iscan in 2:nscans(link)
     meanval += getfield(link.scans[iscan],field)
  end
  return meanval / nscans(link)
end

maxscore1(link::Link) = getmax(link,:score1)
maxscore1(links::Vector{Link}) = maxscore1.(links)

avgscore1(link::Link) = getmean(link,:score1)
avgscore1(links::Vector{Link}) = avgscore1.(links)

maxscore2(link::Link) = getmax(link,:score2)
maxscore2(links::Vector{Link}) = maxscore2.(links)

avgscore2(link::Link) = getmean(link,:score2)
avgscore2(links::Vector{Link}) = avgscore2.(links)

nspecies(link::Link) = link.nspecies
nspecies(links::Vector{Link}) = nspecies.(links)

# Not all data has xic values, so these function have to be special
hasxic(link::Link) = link.hasxic
hasxic(links::Vector{Link}) = hasxic.(links)
maxxic(link::Link) = getmax(link,:xic)
maxxic(links::Vector{Link}) = maxxic.(links)
function avgxic(links::Vector{Link})
  nxic = count(link -> link.hasxic, links)
  meanxic = zeros(nxic)
  ixic = 0
  for link in links
    if link.hasxic
      ixic += 1
      n = 0
      a  = 0.
      for scan in link.scans
        if scan.xic > 0
          n += 1
          a += scan.xic
        end
      end
      meanxic[ixic] = a/n
    end
  end
  return meanxic
end
  
"""

```
count_nspecies(link::Link;mtol=1.0)
```

Count the number of independent species identified in `link`, with a mass difference tolerance
`mtol`, which defaults to `1.0`. This default value is stored in the `link.nspecies` field.

## Example:

```julia-repl
julia> links
 Vector of Links with: 102 links.


julia> count_nspecies(links[1],mtol=1)
7

```

"""
function count_nspecies(link::Link;mtol=1.0)
  nscans(link) == 0 && return 0
  unique = Float64[link.scans[1].mplush]
  for i in 2:nscans(link)
    if findfirst( m -> abs(link.scans[i].mplush-m) < mtol, unique) == nothing
      push!(unique,link.scans[i].mplush)
    end
  end
  return length(unique)
end

"""

```
read_all(xml_file::String=nothing,
         topolink_log = nothing,
         topolink_input = nothing,
         xic_file_name = nothing,
         domain::UnitRange = -1000:1000)
```

Reads the `xml_file` SimXL file, `topolink` log file, `topolink` input file and,
optionally, a file containing `xic` data for each scan, and returns a vector of `Link`
structures containing all data. Optionally, one can set the `domain` to be considered,
consisting of a range of residue numbers. 

## Example:

```julia_repl

julia> using XLStats

julia> links = read_all(xml_file="./data/salbiii_hitsDetail.dat",
                        topolink_log="./data/salbiii_topolink.log",
                        topolink_input="./data/topolink.inp",
                        #xic_file_name="./data/salbiii_xic.dat",
                        domain=2:134)
 Reading all data files... 
 Reading XML file ... 
 Reading XIC data... 
 Vector of Links with: 102 links.

```

"""
function read_all(;xml_file = nothing,
                   topolink_log = nothing,
                   topolink_input = nothing,
                   xic_file_name = nothing,
                   domain::UnitRange = -1000:1000)

  println(" Reading all data files... ")

  error = false
  if xml_file == nothing  
    println("ERROR: You need to provide the xml_file.")
    error = true
  end
  if topolink_log == nothing 
    println("ERROR: You need to provide the topolink log file.")
    error = true
  end
  if topolink_input == nothing
    println("ERROR: You need to provide the topolink input file containing linktype information.")
    error = true
  end

  if error 
    return nothing
  end

  # Read xml file from SIM-XL
  links = read_xml(xml_file,xic_file_name,domain)

  # Read topolink input to get the length of the linkers
  for link in links
    link.dmax = getdmax(link,topolink_input)
  end

  # Read topolink log file to get the euclidean and topological distances
  for link in links 
    link.deuc, link.dtop = readlog(link,topolink_log)
  end

  # Set consistency with 0. tolerance:
  for link in links
    link.consistency = consistency(link,tol=0.)
  end

  return links
end

#
# Function that checks if the link belongs to the chosen domain
#
function indomain(name,domain)
  name == "None" && return false
  r = split(name,'-')
  # Remove wrong assignments for which the two residues are the same
  r[1] == r[2] && return false
  # Check if pair is in domain
  i = parse(Int,r[1][2:end])
  j = parse(Int,r[2][2:end])
  if ( i >= domain[begin] &&
       i <= domain[end] ) &&
     ( j >= domain[begin] &&
       j <= domain[end] )
    return true
  else
    return false
  end
end
indomain(link::Link,domain) = indomain(link.name,domain)

#
# Function that reads the XIC file and fills the xic fields
#
function readxic!(links,xic_file_name)
  
  # Not every xic is provided, -1 means it was not
  for link in links
     link.hasxic = false
     for i in 1:nscans(link)
       link.scans[i].xic = -1
     end
  end

  xic_file_name == nothing && return
  xic_file = open(xic_file_name,"r")
  iread = 0
  xic_in_next_line = false
  for line in eachline(xic_file)
    data = split(line)
    if xic_in_next_line
      for link in links, i in 1:nscans(link)
        if iread == link.scans[i].index
          xic = parse(Float64,data[end])  
          if xic > 0 
            link.hasxic = true
            link.scans[i].xic = xic
          end
          xic_in_next_line = false
          break
        end
      end
    end
    if can_be_parsed_as(Int,data[1])
      iread = parse(Int,data[1])
      xic_in_next_line = true
    end
  end
  close(xic_file)

  nothing
end

function can_be_parsed_as(T::DataType,s::AbstractString)
  try 
    parse(T,s)
    return true
  catch
    return false
  end
end

#
# Function that reads a topolink link log and sets link data
#
function readlog(link,logfile_name)

  file = open(logfile_name,"r") 
  for line in eachline(file)
    comment(line) && continue
    if split(strip(line))[1] == "LINK:"
      residue1_name = oneletter[strip(line[8:12])]
      residue1_index = strip(line[15:19])
      residue2_name = oneletter[strip(line[25:29])]
      residue2_index = strip(line[32:36])
      name = residue1_name*residue1_index*'-'*residue2_name*residue2_index
      if ismatch(name,link.name)
        deuc = parse(Float64,line[43:51])
        dtop = parse(Float64,replace(line[52:61],">"=>""))
        close(file)
        return deuc, dtop
      end
    end
  end
  close(file)
  println(" Warning: Could not find LINK data in TopoLink log for link: $(link.name)")

  return -1., -1.
end

#
# Function that reads the linker length from a topolink input file with linktype
# definitions
#
function getdmax(link,linktype_file)
  file = open(linktype_file,"r") 
  for line in eachline(file)
    comment(line) && continue
    if occursin("linktype",line)
      data = split(line)
      residue1_name = oneletter[data[2]]
      residue2_name = oneletter[data[6]]
      dmax = parse(Float64,data[10])
      name = split(link.name,'-')
      if ( name[1][1:1] == residue1_name && name[2][1:1] == residue2_name ) ||
         ( name[1][1:1] == residue2_name && name[2][1:1] == residue1_name ) 
        close(file)
        return dmax
      end
    end
  end
  close(file)
  return -1.
end

#
# Print link in the topolink format
#
function write(link;tol=0.)
  data = split(link.name,"-")
  res1_name = threeletter[data[1][1:1]]
  res2_name = threeletter[data[2][1:1]]
  res1_index = parse(Int,data[1][2:end])
  res2_index = parse(Int,data[2][2:end])
  cons = consistency(link,tol=tol)
  println(@sprintf("%3s %3i %3s %3i %3i %3.2f %3.2f",
                   res1_name,res1_index,res2_name,res2_index,cons,link.dtop,link.dmax))
end


""" 

```
point_biserial(x::Vector{Bool},y::Vector{<:Real})
```

Compute point-biserial correlation (x assumes 0 or 1 values, y is continuous)
x is boolean (True or False for each group)

"""
function point_biserial(x::AbstractVector{Bool},y::AbstractVector{<:Real})

  ndata = length(x)
  g1 = zeros(Bool,ndata)
  g2 = zeros(Bool,ndata)
  for i in eachindex(x) 
    if x[i] 
      g1[i] = true
    else 
      g2[i] = true
    end
  end
  data1 = y[g1]
  data2 = y[g2]

  n1 = length(data1)
  n2 = length(data2)
  if n1 == 0 || n2 == 0 
    return 0.
  end

  avg1 = mean(data1)
  avg2 = mean(data2)
  sd = std(y)

  pbs = ( ( avg1 - avg2 ) / sd ) * sqrt( n1*n2 / ndata^2 )

  return pbs
end

#
# Function that reads the xml file generated by sim-xl
#
function read_xml(data_file_name,xic_file_name,domain)

  println(" Reading XML file ... ")

  # Count the number of links listed
  data_file = open(data_file_name,"r")
  nlinks = 0
  for line in eachline(data_file)
    comment(line) && continue
    if newlink(line)
      name = setname(line)
      if indomain(name,domain) 
        nlinks = nlinks + 1
      end
    end
  end
  seek(data_file,0)

  # Number of scans per link
  names = String[ "None" for i in 1:nlinks ]
  nscans = zeros(Int,nlinks)

  # Count the number of scans per link
  ilink = 0
  consider_link = false
  for line in eachline(data_file)
    comment(line) && continue
    if newlink(line)
      name = setname(line)
      consider_link = indomain(name,domain)
      if consider_link
        ilink += 1
        names[ilink] = name
      end
    end
    if consider_link && occursin("Scan",line)
      nscans[ilink] += 1
    end
  end
  seek(data_file,0)

  # Initialize the scans per link
  links = [ Link(names[i],nscans[i]) for i in 1:nlinks ]

  # Read the score parameters for each scan
  ilink = 0
  iscan = 0
  for line in eachline(data_file)
    comment(line) && continue
    if newlink(line)
      name = setname(line)
      consider_link = indomain(name,domain)
      if consider_link
        ilink += 1
        iscan = 0
      end
    end
    if consider_link && occursin("Scan",line)
      iscan = iscan + 1

      for pair in data_list
        line = replace(line,pair)
      end
      data = split(line)

      for field in values(data_list)
        ival = findfirst(isequal("$field"),data)
        if ! isnothing(ival) && length(data) > ival
          if field == :peptide
            setfield!(links[ilink].scans[iscan],:pep1,String(data[ival+1]))
            setfield!(links[ilink].scans[iscan],:pep2,String(data[ival+3]))
          else
            if fieldtype(Scan,field) == String
              setfield!(links[ilink].scans[iscan],field,String(data[ival+1]))
            else
              fieldvalue = parse(fieldtype(Scan,field),data[ival+1])
              setfield!(links[ilink].scans[iscan],field,fieldvalue)
            end
          end
        end
      end
    end
  end
  close(data_file)
    
  # Read XIC data
  println(" Reading XIC data... ")
  readxic!(links,xic_file_name)

  # Compute the number of species of each link
  for link in links
    link.nspecies = count_nspecies(link)
  end

  return links
end

#
# Remove a specific link from the list
#
function remove(links,name)
  i=-1
  for link in links
    i=i+1
    if link.name == name
      deleteat!(links,i)
      return links
    end
  end
  error(" Error in remove function: link $(link.name) not found in link list. ") 
end
  

#
# Check if a line is commented
#
comment(line) = (length(strip(line)) == 0) || (strip(line)[1] == '#')

#
# check if line defines a new link
#
newlink(line) = occursin("(",line)
    

#
# Check if a line in a xml file is a new link line, and return
# its name if so
#
function setname(line)
  ! occursin("-",line) && return "None"
  line = replace(line," ("=>"",count=2)
  line = replace(line,") - "=>"-")
  line = replace(line,")"=>"")
  name = strip(line)
  return name
end

function Base.show( io :: IO,::MIME"text/plain", links::Vector{Link} )
  println(" Vector of Links with: ", length(links)," links.")
end

end # module
