module XLStats

using Statistics
using Parameters
using Printf

export Link, getscore, read_all, filter, point_biserial

#
# Convert one letter and three-letter amino acid codes
#
oneletter = Dict("CYS"=>"C", "ASP"=>"D", "SER"=>"S", "GLN"=>"Q", "LYS"=>"K",
                 "ILE"=>"I", "PRO"=>"P", "THR"=>"T", "PHE"=>"F", "ASN"=>"N", 
                 "GLY"=>"G", "HIS"=>"H", "LEU"=>"L", "ARG"=>"R", "TRP"=>"W", 
                 "ALA"=>"A", "VAL"=>"V", "GLU"=>"E", "TYR"=>"Y", "MET"=>"M")
threeletter = Dict( v => k for (k,v) in oneletter )

#
# Structure that contain each link data
#
@with_kw mutable struct Link

  name::String = "None"
  residue1_index::Int = 0
  residue2_index::Int = 0

  consistency::Bool = false
  deuc::Float64 = 0.
  dtop::Float64 = 0.
  dmax::Float64 = 0.

  nscans::Int = 0
  nspecies::Int = 0
  iscan::Vector{Int} = zeros(Int,nscans)
  score1::Vector{Float64} = zeros(nscans)
  score2::Vector{Float64} = zeros(nscans)
  mplush::Vector{Float64} = zeros(nscans)

  hasxic::Bool = false
  xic::Vector{Float64} = zeros(nscans)

end

function Link(name::String,nscans::Int)  
  str = split(name,'-')
  i = parse(Int,str[1][2:end])
  j = parse(Int,str[2][2:end])
  Link(name=name, nscans=nscans, residue1_index=i, residue2_index=j)
end

# 
# syntax sugars for getfield on array of links
#
name(links::Vector{Link}) = getfield.(links,:name)
consistency(links::Vector{Link}) = getfield.(links,:consistency)
deuc(links::Vector{Link}) = getfield.(links,:deuc)
dtop(links::Vector{Link}) = getfield.(links,:dtop)
dmax(links::Vector{Link}) = getfield.(links,:dmax)
nscans(links::Vector{Link}) = getfield.(links,:nscans)
maxscore1(links::Vector{Link}) = maximum.(getfield.(links,:score1))
avgscore1(links::Vector{Link}) = mean.(getfield.(links,:score1))
maxscore2(links::Vector{Link}) = maximum.(getfield.(links,:score2))
avgscore2(links::Vector{Link}) = mean.(getfield.(links,:score2))
nspecies(links::Vector{Link}) = getfield.(links,:nspecies)

# Not all data has xic values, so these function have to be special
maxxic(links::Vector{Link}) = maximum.(getfield.(links,:xic))
function avgxic(links::Vector{Link})
  nxic = count(link -> link.hasxic, links)
  meanxic = zeros(nxic)
  ixic = 0
  for link in links
    if link.hasxic
      ixic += 1
      n = 0
      a  = 0.
      for xic in link.xic
        if xic > 0
          n += 1
          a += xic
        end
      end
      meanxic[ixic] = a/n
    end
  end
  return meanxic
end
  
export name, consistency, deuc, dtop, dmax, nscans, 
       maxscore1, avgscore1,
       maxscore2, avgscore2,
       maxxic, avgxic, 
       nspecies, count_nspecies

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
  @unpack mplush, nscans = link
  unique = Float64[mplush[1]]
  for i in 2:nscans
    if findfirst( m -> abs(mplush[i]-m) < mtol, unique) == nothing
      push!(unique,mplush[i])
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
    link.consistency = setconsistency(link,tol=0.)
  end

  return links
end

#
# Function that compares two links
#
function same_residues(name1,name2)
  r1 = split(name1,'-')
  r2 = split(name2,'-')
  if ( ( r1[1] == r2[1] && r1[2] == r2[2] ) ||  
       ( r1[2] == r2[1] && r1[1] == r2[2] ) )
    return true
  else
    return false
  end
end

#
# Function that checks if the link belongs to the chosen domain
#
function in_domain(name,domain)
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

#
# Function that reads the XIC file and fills the xic fields
#
function readxic!(links,xic_file_name)
  
  # Not every xic is provided, -1 means it was not
  for link in links
     link.xic .= -1
  end

  xic_file_name == nothing && return
  xic_file = open(xic_file_name,"r")
  iread = 0
  xic_in_next_line = false
  for line in eachline(xic_file)
    data = split(line)
    if xic_in_next_line
      for link in links, i in 1:link.nscans
        if iread == link.iscan[i]
          link.xic[i] = parse(Float64,data[end])
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
      deuc = parse(Float64,line[43:51])
      dtop = parse(Float64,replace(line[52:61],">"=>""))
      if same_residues(name,link.name)
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
    if ! comment(line)
      if occursin("linktype",line)
        data = split(line)
        residue1_name = oneletter[data[2]]
        residue2_name = oneletter[data[6]]
        dmax = parse(Float64,data[10])
        name = split(link.name,'-')
        if ( name[1][1] == residue1_name && name[2][1] == residue2_name ) ||
           ( name[1][1] == residue2_name && name[2][1] == residue1_name ) 
          close(file)
          return dmax
        end
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
  consistency = setconsistency(link,tol=tol)
  println(@sprintf("%3s %3i %3s %3i %3i %3.2f %3.2f",
                   res1_name,res1_index,res2_name,res2_index,consistency,link.dtop,link.dmax))
end

#
# Set consistency, given a tolerance
#
function setconsistency(link::Link;tol=0.)
  if (link.dtop >= 0. && link.dtop <= link.dmax + tol)
    return true
  else 
    return false
  end
end

# 
# Compute point-biserial correlation (x assumes 0 or 1 values, y is continuous)
# x is boolean (True or False for each group)
#
function point_biserial(x,y)

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
      if in_domain(name,domain) 
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
      consider_link = in_domain(name,domain)
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
      consider_link = in_domain(name,domain)
      if consider_link
        ilink += 1
        iscan = 0
      end
    end
    if consider_link && occursin("Scan",line)
      iscan = iscan + 1
    
      line = replace(line,"Scan:"=>" ")
      line = replace(line,"Secondary Score:"=>" ")
      line = replace(line,"Score:"=>" ")
      line = replace(line,"Experimental M+H:"=>" ")
      data = split(line)
    
      links[ilink].score1[iscan] =  parse(Float64,(data[2]))
      links[ilink].score2[iscan] =  parse(Float64,(data[3]))
      links[ilink].mplush[iscan] =  parse(Float64,(data[4]))
      links[ilink].iscan[iscan] = parse(Int,data[1])
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
