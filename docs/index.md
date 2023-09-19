VARNA API is a Python interface for [VARNA](http://varna.lri.fr/index.php) (v3-93), a Java lightweight component and applet for drawing the RNA secondary structure.
VARNA allows users to produce drawing in a non-iteractive way via command line.
However, the command line might be massive and complicate in some use cases.
VARNA API aims to simplify such process.
The [online documentation](https://htyao.gitlab.io/varna-api/) is available.

!!! danger "Starting from v1.1.0, VARNA API uses _1-indexed_ to count RNA bases. The change aims to better align with VARNA and ViennaRNA."

## Example

The command below highlights region 11-21 and adds a non-canonical base pair at position (14,20)
on secondary structure `((((((.((((((........)))))).((((((.......))))))..))))))`.
```bash
java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "                                                       " -structureDBN "((((((.((((((........)))))).((((((.......))))))..))))))" -o example.png -auxBPs "(14,20):color=#ff00ff" -highlightRegion "11-21"
```

The equivalence using VARNA API would be
```python
from varnaapi import Structure
v = Structure(structure="((((((.((((((........)))))).((((((.......))))))..))))))")
v.add_highlight_region(11,21)
v.add_aux_BP(14, 20, edge5="s", color="#FF00FF")
v.savefig("example.png")
```
## Installation
```
python3 -m pip install varnaapi
```

## Usage
Here, we show the basic usage of varnaapi.
The first thing after importing `varnaapi` is to setup the location of VARNA to use.

!!! note "By default, the library assumes the VARNA v3-93 in the current folder is used (`./VARNAv3-93.jar`)"

```python
import varnaapi
varnaapi.set_VARNA(path_to_VARNA)
```
Each drawing in VARNA is an object of class inherited from [BasicDraw][varnaapi.BasicDraw]. The standard class to draw from given secondary structure or/and RNA sequence is [Structure][varnaapi.Structure].
```python
ss = "((((((.((((((........)))))).((((((.......))))))..))))))"
v = varnaapi.Structure(structure=ss)
```

One can call the member function `#!python BasicDraw.savefig()` with given file name to save the drawing. The format is either `png` or `svg`, that VARNA will determine from the file name.

```python
v.savefig("my_drawing.png")
```

While using jupyter notebook, one can also set the option `show` to `#!python True` to show the drawing in the notebook

```python
v.savefig("my_drawing.png", show=True)
```
or simply use the member function `#!python BasicDraw.show()`.
```py
v.show()
```

The later one creates a temporary file in `png` format for drawing.

### Configuration

In VARNA, one can chose the structure drawing algorithm, change the style for base pair, or hide backbone in drawing etc.
The full list of parameters with default value can be found [here](config).
In library, the easiest way to modify these parameters is through the member function `#!python BasicDraw.update()`.
Some parameters, such as algorithm, can be set up via specific function. The rest will be supported in future update.

```python
v.update(algorithm='naview', bpStyle='none', drawBackbone=False, bp='#006400')
```

!!! note "Color parameters"
    Value for all color parameters in VARNA API should be readable by the object [colour.Color](https://github.com/vaab/colour), such as human color name, hex, RGB etc.

#### Save configuration

The library offers function to save the parameters setting in YAML format for the further use.
```python
v.dump_param('config.yml')
```

#### Load configuration

There are two ways to load the saved configuration.
The first one is using the member function of `BasicDraw`.
```python
v.load_param('config.yml')
```
The second way is loading the configuration as a global setting in the opened session.
```python
varnaapi.load_config('config.yml')
```

!!! note "Parameter value priority order"
    VARNA API uses parameter values in the following order, _i.e._ if the value in the current order is undefined, then the next one is used

    - Value set by `varnaapi.BasicDraw.update()` or similar functions 
    - Object Value loaded using `varnaapi.BasicDraw.load_param()`
    - Global values loaded using `varnaapi.load_cofig()`
    - Default parameter values.


#### Operations
VARNA allows different operations on drawing, such as highlighting a region, adding auxiliary base pairs.
This can be achieved by using the proper functions. We invite users to read API for more details.
```python
v.add_highlight_region(1, 6, radius=20)
v.add_aux_bp(1, 10, color='red')
```


## Credits
Please kindly cite VARNA [supporting manuscript](https://doi.org/10.1093/bioinformatics/btp250) if you use VARNA API in your research.
Download [bibtex](https://gitlab.inria.fr/amibio/varna-api/-/blob/master/varna.bib) format.
```
Darty, K., Denise, A., & Ponty, Y. (2009). VARNA: Interactive drawing and editing of the RNA secondary structure. Bioinformatics, 25(15), 1974.
```
