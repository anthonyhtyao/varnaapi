VARNA API is a Python interface for [VARNA](http://varna.lri.fr/index.php), a Java lightweight component and applet for drawing the RNA secondary structure.
VARNA allows users to produce drawing in a non-iteractive way via command line.
However, the command line might be massive and complicate in some use cases.
VARNA API aims to simplify such process.

__NOTE__: The VARNA API is 0-indexed unlike VARNA, which is 1-indexed.

## Example

The command below highlights region 11-21 and adds a non-canonical base pair at position (14,20)
on secondary structure `((((((.((((((........)))))).((((((.......))))))..))))))`.
```bash
java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "                                                       " -structureDBN "((((((.((((((........)))))).((((((.......))))))..))))))" -o example.png -algorithm radiate -auxBPs "(14,20):color=#FF00FF,thickness=1.0,edge5=s,edge3=wc,stericity=cis" -highlightRegion "11-21:radius=15.0,fill=#9999FF,outline=#3333FF"
```

The equivalence in python is
```python
from varnaapi import VARNA
v = VARNA(structure="((((((.((((((........)))))).((((((.......))))))..))))))")
v.add_highlight_region(10,20)
v.add_aux_BP(13, 19, edge5="s", color="#FF00FF")
v.savefig("example.png")
```
## Installation
```
python3 -m pip install varnaapi
```

## Usage
Here, we show the basic usage of varnaapi. Please refer the [API](https://htyao.gitlab.io/varna-api/varna) page for more details.
The first thing after importing `varnaapi` is to setup the location of VARNA used.
The default is `VARNAv3-93.jar` in the current folder.
```python
import varnaapi
varnaapi.set_VARNA(path_to_VARNA)
```
Each drawing in VARNA is an object called `VARNA` created from given secondary structure or/and RNA sequence.
```python
ss = "((((((.((((((........)))))).((((((.......))))))..))))))"
v = varnaapi.VARNA(structure=ss)
```
Then we can add operations on drawing by calling object functions, such as `VARNA.set_algorithm()` to choose a drawing algorithm, `VARNA.add_highlight_region()` to highlight a region etc. 
```python
v.set_algorithm('line')
v.add_highlight_region(0, 5, radius=20)
```
Finally, we can draw the secondary structure
```python
v.savefig(path_to_store)
```

## Credits
Please kindly cite VARNA [supporting manuscript](https://doi.org/10.1093/bioinformatics/btp250) if you use VARNA API in your research.
Download [bibtex](https://gitlab.inria.fr/amibio/varna-api/-/blob/master/varna.bib) format.
```
Darty, K., Denise, A., & Ponty, Y. (2009). VARNA: Interactive drawing and editing of the RNA secondary structure. Bioinformatics, 25(15), 1974.
```
