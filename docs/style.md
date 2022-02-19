# Style classes for VARNA operations

Two helper class groups can be imported from `#! python varnaapi.param` for two VARNA operations, [add_annotation][varnaapi.BasicDraw.add_annotation] and [add_bases_style][varnaapi.BasicDraw.add_bases_style].

## Annotation

An `Annotation` object represents a textual annotation added to a VARNA drawing.
The object stores the text and other informtation needed.
One can add `Annotation` to drawing using [BasicDraw.add_annotation][varnaapi.BasicDraw.add_annotation].
Four annotation types allowed in VARNA are represented by four objects below.


::: varnaapi.param
    selection:
		  members: ["BaseAnnotation", "LoopAnnotation", "HelixAnnotation", "StaticAnnotation"]
    rendering:
      heading_level: 3

::: varnaapi.param
    selection:
		  members: ["BasesStyle"]

