# XLStats

A package to analyze the correlations between structural distances and mass-spectrometry quantitative data of chemical-crosslinks obtained with SimXL. The structural properties are supposed to be computed from models (PDB files) using TopoLink.

## Installing

```julia
julia> ] add https://github.com/m3g/XLStats.jl
```

## Example

The example files are available at the `test/data` directory. 

### Reading the data

```julia
julia> using XLStats

julia> links = read_all(xml_file="./data/salbiii_hitsDetail.dat",
                        topolink_log="./data/salbiii_topolink.log",
                        topolink_input="./data/topolink.inp",
                        xic_file_name="./data/salbiii_xic.dat",
                        domain=2:134)
 Reading all data files... 
 Reading XML file ... 
 Reading XIC data... 
 Vector of Links with: 102 links.
```

### Available parameters

The goal is to probe correlations between the spectral parameters and the topological (surface-accessible) distances obtained with TopoLink. The available parameters to be compared are (fields of the `Link` structure, which are the elements of the `links` vector created above):

| field name    | Meaning                       | type of value     | Example                   |
|:-------------:|:-----------------------------:|:-----------------:|:-------------------------:|
| `name1`       | Name of link                  | `String`          | `"K17-K6"`                |
| `consistency` | Consitency with the structure | `Bool`            | `false`                   |
| `deuc`        | Euclidean distance            | `Float64`         | `16.696`                  |
| `dtop`        | Toplogical distance           | `Float64`         | `18.218`                  |
| `dmax`        | Maximum linker length         | `Float64`         | `17.000`                  |
| `nscans`      | Number of scans of the XL     | `Int64`           | `15`                      |
| `nspecies`    | Number of species of the XL   | `Int64`           | `7`                       |
| `score1`      | `score1` of the scans         | `Vector{Float64}` | `[6.104, 5.978,... ]`     |
| `score2`      | `score2` of the scans         | `Vector{Float64}` | `[3.379, 2.668,... ]`     |
| `hasxic`      |  XIC data is available?       | `Bool`            | `false`                   |
| `xic`         | XIC data for each scan        | `Vector{Float64}` | `[-1.0, 157034.76,... ]`  |

### Ploting correlations

For example, lets plot the correlation between the topological distance of each link and the maximum `score1`. There are helper functions (see below) for each of these scores:

```julia
julia> using Plots

julia> x = deuc(links);

julia> y = maxscore1(links);

julia> scatter(x,y,label="",xlabel="Euclidian distance",ylabel="Maximum score1")

```
produces:

![fig1](./test/plots/score1_vs_deuc.png)

Note that the functions `deuc` and `maxscore1`, in the example, return vectors of data extracted for every link in the `links` list. Similar functions exist for the other parameters, and are:

`name`, `resnames`, `indexes, `indomain`,
`consistency`, `deuc`, `dtop`, `dmax`, `nscans`,
`maxscore1`, `avgscore1`, `maxscore2`, `avgscore2`, `maxxic`, `avgxic`, `nspecies`, `count_nspecies`

*Note:* The `consistency` function is special. It returns `true` or `false` depending if the topological distance (`dtop`) of the XL is smaller than the maximum linker reach (`dmax`), up to a tolerance `tol`. By default `tol=0`, but the function accepts the tolerance as an argument. Thus, for example:  

```julia
julia> x = consistency(links)
```
will return all links for which `dtop < dmax`. But  
```julia
julia> x = consistency(links,tol=2.0)
```
will return all links for which `dtop < dmax + 2.0`.

### Filter links

The functions above allow filtering the links by any criteria, using the standard `filter` function. For example,

```julia
julia> bad_links = filter( link -> consistency(link,tol=2.0) == false, links )
 Vector of Links with: 80 links.

julia> large_deuc = filter( link -> deu(link) > 10, links )
 Vector of Links with: 90 links.

julia> large_score1 = filter( link -> maxscore1(link) > 5., links )
 Vector of Links with: 11 links.

```

The `resnames` and `indexes` functions also allow filtering by residue types and indexes. 
Since the permutation between types are usually meaningless, the `ismatch` function is provided:

```julia
julia> resnames(links[2])
("SER", "LYS")

julia> ismatch(resnames(links[2]),("SER","LYS"))
true

julia> ismatch(resnames(links[2]),("LYS","SER"))
true

julia> indexes(links[2])
(131, 99)

julia> ismatch(indexes(links[2]),(99,131))
true
```
Thus this function can be used for proper filtering of links by type, for example:

```julia
julia> filter( link -> ismatch(resnames(link),("SER","LYS")), links )
 Vector of Links with: 20 links.
```

To filter the links by residue numbers, the `indomain` function is mostly useful:
```julia
julia> indexes(links[2])
(131, 99)

julia> indomain(links[2],1:131)
true

julia> indomain(links[2],1:120)
false

julia> mydomain = filter(link -> indomain(link,1:120), links)
 Vector of Links with: 77 links.
```

### Remarks on XIC data

XIC data might not be available for every XL, or every scan of every XL. Scans without XIC data will show `-1` at the XIC field. The links for which XIC data is available can be filtered with:

```julia
julia> links_with_xic = filter(link -> link.hasxic, links)
 Vector of Links with: 32 links.

```

and then it makes sense to use the `avgxic` and `maxxic` functions in this subset of data. Note that `avgxic` will compute the average value of XIC data disregarding any missing XIC, even if not every scan for the link has XIC data. 

XIC data also spans many orders of magnitude. Thus, converting it a `log10` scale is handy. This can be readily using the `.` (dot) syntax in Julia, for example:

```julia

julia> x = dtop(links_with_xic);

julia> y = log10.(maxxic(links_with_xic));

julia> scatter(x,y,label="",xlabel="Topological distance",ylabel="log₁₀(maxxic)")

```

(note the `.` after `log10`, which means that the function will be applied to every element of the array that follows). 
This produces:

![fig2](./test/plots/logmxic_vs_top.png)

### Exporting data 

To export data to be analyzed by other software use the `DelimitedFiles` package:

- Install with
```julia
julia> ] add DelimitedFiles
```

- For example, write a data file with
```julia
julia> using DelimitedFiles

julia> writedlm("mydata.txt",[x y])
```

where `x` and `y` could be obtained by the commands of the examples above.

















