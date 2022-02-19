In VARNA API, all drawing classes are inherited from the class [BasicDraw][varnaapi.BasicDraw]. We offer three Python classes [Structure][varnaapi.Structure], [Comparison][varnaapi.Comparison], and [Motif][varnaapi.Motif] for different RNA secondary structure visualisation.
The first two correspond to the classic and the comparison mode in VARNA. The last one is the special case for motif drawing.

- [Structure][varnaapi.Structure]: The most common drawing class taking a secondary structure and/or an RNA sequence as input.
- [Comparison][varnaapi.Comparison]: Comparison between two structures and two sequences.
- [Motif][varnaapi.Motif]: Customized class to draw motif, an union of loops.

::: varnaapi.models
    selection:
      members: ["BasicDraw"]
      inherited_members: True

::: varnaapi.models
    selection:
      members: ["Structure"]

::: varnaapi.models
    selection:
      members: ["Comparison"]

::: varnaapi.models
    selection:
      filters: ["!savefig"]
      members: ["Motif"]

![Motif ((*)(*)(((*)(*))))](assets/images/motif_ex.png)

<figcaption>Motif ((*)(*)(((*)(*))))</figcaption>

### Example
Figure above is created with
```python
from varnaapi import Motif
from varnaapi.param import BaseAnnotation
m = Motif("((*)(*)(((*)(*))))", sequence="  *AU* CC *  *    ")
m.add_annotation(BaseAnnotation(" Root", 1))
m.add_annotation(BaseAnnotation("Dummy", 13))
# Show how base indices work for motif.
# Remeber that VARNA is 1-indexed
m.set_default_color(baseNum="#a6a6a6")
m.set_numeric_params(periodNum=4)
m.savefig("motif_ex.png")
```
