# XLStats

A package to analyze the correlations between structural distances and mass-spectrometry quantitative data of chemical-crosslinks obtained with SimXL. The structural properties are supposed to be computed from models (PDB files) using TopoLink.

## Installing

```
julia> ] add https://github.com/m3g/XLStats.jl

```

## Example

The example files are available at the `test/data` directory. 

```julia-repl

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

