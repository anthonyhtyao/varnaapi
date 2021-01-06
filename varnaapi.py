import re
import os
import abc
from string import ascii_lowercase, ascii_uppercase

__version__ = '0.1.0'

_VARNA_PATH="VARNAv3-93.jar"
HEX = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')
BORDER = re.compile('^\d+x\d+$')

DEFAULT_COLOR_LIST = ['backbone', 'background', 'baseInner', 'baseName', 'baseNum',
                      'baseOutline', 'bp', 'gapsColor', 'nsBasesColor']
"""Allowed options for color setting

| Name         | Object in panel                                         |
|--------------|---------------------------------------------------------|
| backbone     | Phosphate-sugar backbone (aka skeleton) of the RNA      |
| background   | Background color used within the panel                  |
| baseInner    | Inner base color                                        |
| baseName     | Nucleotide name                                         |
| baseNum      | Base numbers                                            |
| baseOutline  | Outer base color                                        |
| bp           | Base-pair                                               |
| nsBasesColor | Non-standard bases (Anything but `A`, `C`, `G` or `U`)  |
"""

OPTIONS = ['autoHelices', 'autoInteriorLoops', 'autoTerminalLoops', 'drawBackbone', 'drawBases', 'drawNC', 'drawTertiary', 'fillBases', 'flat']
"""Boolean option list

| Name              | Option                                                                                                                                        | Default |
|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|---------|
| autoHelices       | Annotates each and every helix in the RNA with a unique `Hn` label                                                                            | False   |
| autoInteriorLoops | Annotates each and every interior loop of the RNA with a unique `In` label                                                                    | False   |
| autoTerminalLoops | Annotates each and every terminal loop of the RNA with a unique `Tn` label                                                                    | False   |
| drawBackbone      | Backbone drawing                                                                                                                              | True    |
| drawBases         | Displays the outline of a nucleotide base                                                                                                     | True    |
| drawNC            | Displays non-canonical base-pairs                                                                                                             | True    |
| drawTertiary      | Display of `non-planar` base-pairs, _i.e._ pseudoknots [^1]                                                                                   | False   |
| fillBases         | Fill bases                                                                                                                                    | True    |
| flat              | In `radiate` drawing mode, redraws the whole structure, aligning to a  baseline the base featured on the exterior loops (aka "dangling ends") | False   |

[^1]: Since there is no canonical definition of pseudoknotted portions, a maximal planar subset is extracted from the input structure, defined to be the planar portion, and used as a scaffold for the drawing algorithms.

"""

NUMERIC_PARAMS = ['bpIncrement', 'periodNum', 'resolution', 'rotation', 'spaceBetweenBases', 'zoom']
"""Allowed numeric parameters

| Label             | Type  | Description                                                                                                                                                                                                    |
|-------------------|-------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| bpIncrement       | float | In linear drawing mode, defines the vertical increment used to separate two successive, nested base-pairs                                                                                                      |
| periodNum         | int   | Sets the interval between two successive base numbers. More specifically, if `k` is the period, then the first and last bases  of the RNA are numbered, along with each base whose number is a multiple of `k` |
| resolution        | float | Chooses the resolution of a bitmap export, _i.e._ the multiplier in  the number of pixels in each dimension of the exported picture.                                                                           |
| rotation          | float | Rotates the whole RNA of a certain angular increment                                                                                                                                                           |
| spaceBetweenBases | float | Sets the distance between consecutive bases                                                                                                                                                                    |
| zoom              | float | Defines the level of zoom and zoom increment used to display the RNA within this panel                                                                                                                         |
"""

BP_STYLES = ['none', 'simple', 'rnaviz', 'lw']
"""Allowed options for base-pair style

| Label  | Description                                                                                            |
|--------|--------------------------------------------------------------------------------------------------------|
| none   | Base-pairs are not drawn, but can be implicitly seen from "ladders", _i.e_ helix structures            |
| simple | A simple line is used to draw any base-pair, regardless of its type                                    |
| rnaviz | A small square is drawn at equal distance of the two partners                                          |
| lw     | Both canonical and non-canonical base-pairs are rendered according to the Leontis/Westhof nomenclature |

__See Also:__ [VARNA.set_bp_style][varnaapi.VARNA.set_bp_style]

"""

PARENTHESES_SYSTEMS = [
    ("(", ")"),
    ("[", "]"),
    ("<", ">"),
    ("{", "}")
] + [(c1, c2) for c1, c2 in zip(ascii_uppercase, ascii_lowercase)]
PARENTHESES_OPENING = [c1 for c1, c2 in PARENTHESES_SYSTEMS]
PARENTHESES_CLOSING = {c2: c1 for c1, c2 in PARENTHESES_SYSTEMS}

def set_VARNA(path):
    """Set VARNA location
    """
    global _VARNA_PATH
    _VARNA_PATH = path

class BasesStyle:
    """Defines a custom base-style, to be applied later to a set of bases.
    A BasesStyle style contains colors used for different components of a base.
    See [\_\_init\_\_][varnaapi.BasesStyle.__init__] for more details.

    __See Also:__ [VARNA.add_bases_style][varnaapi.VARNA.add_bases_style]
    """
    def __init__(self, fill=None, outline=None, label=None, number=None):
        """Basesstyle constructor from given colors for different components.
        At least one argument should be given.

        Args:
            fill (Hex): color of inner part of base
            outline (Hex): color of outline of base
            label (Hex): base text color
            number (Hex): base number color

        Examples:
            >>> style = BasesStyle(fill='#FF0000', outline='#00FF00')
        """
        self.color = {}
        self.update(fill, outline, label, number)

    def update(self, fill=None, outline=None, label=None, number=None):
        """Update component colors.
        Same rule as [\_\_init\_\_][varnaapi.BasesStyle.__init__]
        """
        if fill is None and outline is None and label is None and number is None:
            raise Exception("At least one should not be None")
        if fill is not None:
            fill = assert_hex_color(fill)
            self.color["fill"] = fill
        if outline is not None:
            outline = assert_hex_color(outline)
            self.color["outline"] = outline
        if label is not None:
            label = assert_hex_color(label)
            self.color["label"] = label
        if number is not None:
            number = assert_hex_color(number)
            self.color["number"] = number

    def __str__(self):
        order = ['fill', 'outline', 'label', 'number']
        lst = ["{}={}".format(k, self.color[k]) for k in order if k in self.color]
        return ",".join(lst)


class _Annotation:
    """Basic Annotation
    """
    def __init__(self, text, type, color="#000000", size=10):
        self.text = text
        self.type = type
        self.color = color
        self.size = size #: int: font size

    def asdict(self):
        return {'text': self.text, 'type': self.type, 'color': self.color,
                'size': self.size}

    @abc.abstractmethod
    def to_cmd(self):
        pass


class _ObjectAnnotation(_Annotation):
    def __init__(self, text, type, anchor, color="#000000", size=10):
        super().__init__(text, type, color, size)
        self.anchor = anchor

    def asdict(self):
        d = super().asdict()
        d['anchor'] = self.anchor + 1
        return d

    def to_cmd(self):
        return "{text}:type={type},anchor={anchor},color={color},size={size}"\
            .format(**self.asdict())


class BaseAnnotation(_ObjectAnnotation):
    def __init__(self, text:str, anchor:int, color="#000000", size=10):
        """Annoation on a base.

        Args:
            text: Annotation caption
            anchor: Index of base to annotate
            color (Hex): Annotation color
            size (int): Font size
        """
        super().__init__(text, 'B', anchor, color, size)


class LoopAnnotation(_ObjectAnnotation):
    """Same as [BaseAnnotation][varnaapi.BaseAnnotation] but on a loop.
    Argument `anchor` can be index of any base in the loop of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=10):
        super().__init__(text, 'L', anchor, color, size)


class HelixAnnotation(_ObjectAnnotation):
    """Same as [BaseAnnotation][varnaapi.BaseAnnotation] but on an helix.
    Argument `anchor` can be index of any base in the helix of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=10):
        super().__init__(text, 'H', anchor, color, size)


class StaticAnnotation(_Annotation):
    def __init__(self, text, x, y, color="#000000", size=10):
        """Annotation on a specified position in VARNA drawing.
        Unlike [BaseAnnotation][varnaapi.BaseAnnotation], argument `anchor` is omitted.
        However, arguments `x` and `y` are needed to specify annotation position.

        __Note:__ It is unrecommended to use static annotation unless you know what you're doing

        Args:
            x (int): x-coordinate of position
            y (int): y-ccordinate of position

        Examples:
            >>> sa = StaticAnnotation("Hello World", 100, 150, color="#FF0000")
        """
        super().__init__(text, 'P', color, size)
        self.x = x
        self.y = y

    def asdict(self):
        d = super().asdict()
        d['x'] = self.x
        d['y'] = self.y
        return d

    def to_cmd(self):
        return "{text}:type={type},x={x},y={y},color={color},size={size}"\
            .format(**self.asdict())

def is_hex_color(color):
    match = HEX.search(color)
    if match:
        return True
    return False

def assert_hex_color(color):
    if not is_hex_color(color):
        raise Exception("{} is not a valid Hex color code".format(color))
    return color.upper()

def assert_is_number(name, value):
    if not (isinstance(value, float) or isinstance(value, int)):
        raise Exception(name + " should be a float or an integer")

def assert_valid_interval(length, *args):
    for i in args:
        if i < 0 or i >= length:
            raise Exception("{} out of range".format(args))

def check_structure(ss):
    pass


def _bp_to_struct(bps):
    """Base pair list to structure"""
    n = max([j for i, j in bps]) + 1
    ss = [-1 for i in range(n)]
    for i, j in bps:
        ss[i], ss[j] = j, i
    return ss

def _parse_vienna(ss):
    """
    Parse secondary structure in dot-bracket notation
    """
    stacks = {c:[] for c in PARENTHESES_OPENING}
    res = [-1 for i in range(len(ss))]
    for i,c in enumerate(ss):
        if c in PARENTHESES_OPENING:
            stacks[c].append(i)
        elif c in PARENTHESES_CLOSING:
            ii = stacks[PARENTHESES_CLOSING[c]].pop()
            res[ii],res[i] = i,ii
    return res

class VARNA:
    def __init__(self, seq:str=None, structure=None):
        """Classic VARNA drawing mode. Constructor from given RNA sequence or/and secondary structure.
        If sequence and structure have different size, the larger one is used
        and ` `s or `.`s will be added to sequence or structure to complement.

        Args:
            seq: Raw nucleotide sequence for the displayed RNA.
                 Each base must be encoded in a single character.
                 Letters others than `A`, `C`, `G`, `U` and space are tolerated.
            structure (str or list): RNA (pseudoknotted) secondary structure in one of three formats

              - Dot-Bracket Notation (DBN)
              - List of pair of int representing a list of base-pairs
              - List of int, in which i-th value is `j` if `(i,j)` is a base pair or `-1` if i-th base is unpaired

        """
        self.length = -1
        self.structure = []
        self.sequence = ""
        self._init_features()

        if structure is not None:
            if isinstance(structure, list):
                if len(structure) > 0:
                    first = structure[0]
                    if len(first)==1:
                        self.structure = check_structure(structure)
                    elif len(first)==2:
                        self.structure = _bp_to_struct(structure)
                    else:
                        raise Exception("Unrecognized structure format for %s"%(structure))
            # Dot-Bracket Notation
            elif isinstance(structure, str):
                self.structure = _parse_vienna(structure)
                self.dbn = structure
            self.length = len(self.structure)
        if seq is not None:
            self.length = max(self.length,len(seq))
            self.sequence = seq
        # Now we know the length, let's extend the sequence and structure if necessary
        self.sequence += " "*(self.length-len(self.sequence))
        self.structure += [-1]*(self.length-len(self.structure))

    def _init_features(self):
        self.aux_BPs = []
        self.highlight_regions = []
        self.params = {'algorithm': "radiate"}
        self.default_color = {}
        self.options = {}
        self.title = None
        self.bases_styles = {}
        self.annotations = []

    def add_aux_BP(self, i:int, j:int, edge5:str='wc', edge3:str='wc', stericity:str='cis', color='#0000FF', thickness:float=1.0):
        """Add an additional base pair `(i,j)`, possibly defining and using custom style

        Args:
            i: 5' position of base pair
            j: 3' position of base pair
            edge5: Edge 5' used for interaction in non-canonical base-pairs, as defined by the Leontis/Westhof classification of base-pairs. Admissible values are __wc__ (Watson/Crick edge), __h__ (Hoogsteen edge) and __s__ (Sugar edge).
            edge3: Edge 3' used for interaction in non-canonical base-pairs. Admissible values are __wc__, __h__ and __s__.
            stericity: Orientation of the strands. Admissible values are __cis__ and __trans__
            color (Hex): Base-pair color in Hex color codes
            thickness: Base-pair thickness
        """
        assert_valid_interval(self.length, i, j)
        if edge5 not in ['wc', 's', 'h']:
            raise Exception("edge5 should be one of wc, s, and h")
        if edge3 not in ['wc', 's', 'h']:
            raise Exception("edge3 should be one of wc, s, and h")
        if stericity not in ['cis', 'trans']:
            raise Exception("stericity should be either cis or trans")
        color = assert_hex_color(color)
        assert_is_number('thickness', thickness)

        self.aux_BPs.append((i+1, j+1, color, thickness, edge5, edge3, stericity))

    def add_highlight_region(self, i:int, j:int, radius:float=15.0, fill="#9999FF", outline="#3333FF"):
        """Highlights a region by drawing a polygon of predefined radius,
        fill color and outline color around it.
        A region consists in an interval from base `i` to base `j`.

        Args:
            i: 5'-end of the highlight
            j: 3'-end of the highlight
            radius: Thickness of the highlight
            fill (Hex): The color used to fill the highlight
            outline (Hex): The color used to draw the line around the highlight
        """
        assert_valid_interval(self.length, i, j)
        fill = assert_hex_color(fill)
        outline = assert_hex_color(outline)
        assert_is_number('radius', radius)

        self.highlight_regions.append((i+1, j+1, radius, fill, outline))

    def set_algorithm(self, algo):
        """Set algorithm other than __naview__ to draw secondary structure.
        Supported options are __line__, __circular__, __radiate__ and __naview__.
        """
        if algo not in ['line', 'circular', 'radiate', 'naview']:
            raise Exception("Sould be one of line, circular, radiate or naview")
        self.params['algorithm'] = algo


    def set_title(self, title:str, color:str='#808080', size:int=10):
        """Set title displayed at the bottom of the panel with color and font size
        """
        self.title = (title, color, size)

    def set_zoom_level(self, level:float):
        """Defines the level of zoom and zoom increment used to display the RNA within this panel"""
        self.params['zoom'] = level

    def set_default_color(self, **kwargs):
        """Set default color used for different objects in the panel.

        Args:
            **kwargs (dict): See [DEFAULT_COLOR_LIST][varnaapi.DEFAULT_COLOR_LIST] for the list of allowed keywords.
                Value of a keyword is the default color, in Hex color codes, used for the related object.

        Examples:
            >>> set_default_color({'backbone': '#000000', 'bp':'#FFFF00'})

        """
        for key, value in kwargs.items():
            if key not in DEFAULT_COLOR_LIST:
                raise Exception("{} is not a valid keyword".format(key))
            assert_hex_color(value)
        self.default_color = kwargs

    def toggle_options(self, **kwargs):
        """Enable or disable options of drawing

        Args:
            **kwargs (dict): See [OPTIONS][varnaapi.OPTIONS] for detailed option lists.
                      Value of keyword is either `True` or `False`


        """
        for key, value in kwargs.items():
            if key not in OPTIONS:
                raise Exception("{} is not a valid keyword".format(key))
            if not isinstance(value, bool):
                raise Exception(key + "should be a boolean")
        self.options = kwargs

    def set_numeric_params(self, **kwargs):
        """Change value of numeric parameters in one function.
        This is equivalent to use setting function of each parameter,
        such as [set_bp_increment][varnaapi.VARNA.set_bp_increment].

        Args:
            **kwargs (dict): See [NUMERIC_PARAMS][varnaapi.NUMERIC_PARAMS] for allowed keywords.

        """
        for key, value in kwargs.items():
            if key not in NUMERIC_PARAMS:
                raise Exception("{} is not a valid keyword".format(key))
            assert_is_number(key, value)
            if key == "periodNum":
                self.params['periodNum'] = int(value)
            else:
                self.params[key] = float(value)

    def format_structure(self):
        """Return secondary structure in dot-brackaet notation
        """
        def greedy_fill(c1, c2, res, ss, i, j):
            if i <= j:
                k = ss[i]
                if k == -1:
                    greedy_fill(c1, c2, res, ss, i+1, j)
                elif k > i:
                    if k <= j:
                        res[i], res[k] = c1, c2
                        ss[i], ss[k] = -1, -1
                        greedy_fill(c1, c2, res, ss, i+1, k-1)
                        greedy_fill(c1, c2, res, ss, k+1, j)

        res = ["." for _ in range(self.length)]
        ss = self.structure[:]
        for c1, c2 in PARENTHESES_SYSTEMS:
            greedy_fill(c1, c2, res, ss, i=0, j=self.length-1)
            finished = True
            for i in ss:
                if i != -1:
                    finished = False
            if finished:
                break
        return "".join(res)


    def set_border(self, border:str):
        """Sets the width and height of the panel border, _i.e._ the gap
        between the panel boundaries and those of the surface used to draw the RNA.
        Parameter `border` is in format `"wxh"` where `w` and `h` are width and height separated by `x`.

        Example:
            >>> set_border("20x30")
        """
        match = BORDER.search(border)
        if not match:
            raise Exception("border should be the format nxm where n and m are numbers")
        self.params['border'] = "\"{}\"".format(border)

    def set_bp_style(self, style:str):
        """Set default style for base-pairs rendering, chosen among [BP_STYLES][varnaapi.BP_STYLES]

        __Note:__ `lw` is set by default

        Example:
            >>> varna.set_bp_style("simple")
        """
        if style not in BP_STYLES:
            raise Exception('Should be one of {}'.format(BP_STYLES))
        self.params['bpStyle'] = style

    def set_bp_increment(self, value:float):
        """In linear drawing mode, defines the vertical increment used to
        separate two successive, nested base-pairs.

        Example:
            >>> varna.set_bp_increment(1.2)
        """
        assert_is_number('value', value)
        self.params['bpIncrement'] = float(value)


    def add_bases_style(self, style:BasesStyle, bases:list):
        """Apply a [BasesStyle][varnaapi.BasesStyle] to a list of positions.
        If a position is assigned to more than one styles,
        one of them will be randomly used.

        Args:
            style: Style to apply
            bases: List of 0-indexed positions

        Examples:
            >>> style1 = BasesStyle(fill="#FF0000")
            >>> style2 = BasesStyle(fill="#FFFF00" outline="#00FF00")
            >>> varna.add_bases_style(style1, [0,2,4])
            >>> varna.add_bases_style(setye1, [10,11,12])
            >>> varna.add_bases_style(style2, [4,5,6,7])

        """
        if not isinstance(style, BasesStyle):
            raise Exception("style should be BasesStyle object")
        if len(bases) > 0:
            self.bases_styles[style] = self.bases_styles.get(style, set()).union({i+1 for i in bases})


    def add_annotation(self, annotation:_Annotation):
        """Add an annotation.
        Argument should be a valid [Annotation](annotation.md) object

        Examples:
            >>> a = LoopAnnotation("L1", 6, color="#FF00FF")
            >>> varna.add_annotation(a)
        """
        # Assert is annotation
        if not isinstance(annotation, _Annotation):
            raise Exception("Should be a valid annotation object")
        self.annotations.append(annotation)

    def _gen_command(self):
        """
        Return command to run VARNA
        """
        cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd".format(_VARNA_PATH)

        cmd += self._gen_input_cmd()

        cmd += " -o {}".format(self.output)

        # Command for defualt colors
        for key, color in self.default_color.items():
            if color is not None:
                cmd += " -{} \"{}\"".format(key, color)

        # Options
        for key, value in self.options.items():
            if value is not None:
                cmd += " -{} {}".format(key, value)

        # Params
        for key, value in self.params.items():
            cmd += " -{} {}".format(key, value)

        # Title
        if self.title is not None:
            cmd += " -title {} -titleColor {} -titleSize {}".format(*self.title)

        # Aux Base pairs
        if len(self.aux_BPs) > 0:
            auxbps = ["({},{}):color={},thickness={},edge5={},edge3={},stericity={}"
                      .format(*t) for t in self.aux_BPs]
            cmd += " -auxBPs \"{}\"".format(";".join(auxbps))

        # Highlight Region
        if len(self.highlight_regions) > 0:
            regions = ["{}-{}:radius={},fill={},outline={}".format(*t) for t in self.highlight_regions]
            cmd += " -highlightRegion \"{}\"".format(";".join(regions))

        # BasesStyles
        for ind, (style, bases) in enumerate(self.bases_styles.items()):
            cmd += " -basesStyle{} {}".format(ind+1, str(style))
            cmd += " -applyBasesStyle{}on {}".format(ind+1, ','.join(map(str, bases)))

        # Annotations
        if len(self.annotations) > 0:
            cmd += " -annotations \"{}\"".format(';'.join([t.to_cmd() for t in self.annotations]))

        return cmd


    def _gen_input_cmd(self):
        return " -sequenceDBN \"{}\" -structureDBN \"{}\"".format(self.sequence, self.format_structure())

    def savefig(self, output):
        """
        Call VARNA to draw and store the paint in output
        """
        self.output = output
        cmd = self._gen_command()
        print(cmd)
        os.popen(cmd).close()

    def __repr__(self):
        return repr((self.format_structure(),self.sequence))


class Comparison(VARNA):
    def __init__(self, seq1, structure1, seq2, structure2):
        """Drawing of two aligned RNAs.
        Unlike classic [VARNA][varnaapi.VARNA] mode,
        both sequences and structures __MUST__ be specified and have the same size.
        Additionally, the merged secondary structures must currently be without any crossing
        interactions (e.g. pseudoknots), and merging them should give a secondary structure.
        Gap character is `.`.
        Args:
            seq1 (str): Sets the gapped nucleotide sequence for the first RNA sequence
            structure1 (str): Sets the first secondary structure in Dot-Bracket Notation
            seq2 (str): Sets the gapped nucleotide sequence for the second sequence
            strcuture2 (str): Sets the second secondary structure in Doc-Bracket Notation
        """
        if not (len(seq1) == len(structure1) == len(seq2) == len(structure2)):
            raise Exception("All length should be equal")
        self.seq1 = seq1
        self.structure1 = structure1
        self.seq2 = seq2
        self.structure2 = structure2
        self.length = len(seq1)
        self._init_features()

    def _gen_input_cmd(self):
        return " -comparisonMode True -firstSequence \"{}\" -firstStructure \"{}\" -secondSequence \"{}\" -secondStructure \"{}\"".format(self.seq1, self.structure1, self.seq2, self.structure2)

    def __repr__(self):
        return repr((self.seq1, self.structure1, self.seq2, self.structure2))


class Motif(VARNA):
    def __init__(self, motif, sequence=None):
        """Special class for motif drawing.
        A motif is a rooted ordered tree, similar to a secondary structure,
        but whose leaves may represent base paired positions, named open base
        pairs or open paired leaves and denoted by `(*)`, and the root always
        represents a closing base pair. A motif can also be seen as an union
        of consecutive loops. The figure below represents `((*)(*)(((*)(*))))`.

        Motif class inherits from [VARNA][varnaapi.VARNA] with some pre-set
        parameters.

        - rotation is set at `180`
        - default base pair style is `simple`
        - base number is hidden by setting default color to white
        (default background color)

        A dummy base pair is added after each open base pair and in front of
        the root, as shown in the figure below.
        Therefore, the index of bases is changed after creating the object.
        For example, the index of first base of root is `1` instead of `0`.
        The default bases style for root is
        `BasesStyle(fill="#606060", outline="#FFFFFF",number="#FFFFFF")` and
        `BasesStyle(fill="#DDDDDD", outline="#FFFFFF", number="#FFFFFF")` for
        dummy bases. One can change them using
        [set_root_bases_style][varnaapi.Motif.set_root_bases_style] and
        [set_dummy_bases_style][varnaapi.Motif.set_dummy_bases_style].

        Args:
            motif (str): Motif in Dot-Bracket Notation.
                `(*)` is used to represent open base pair.
            sequence (str): Chain of characters for motif. Note that sequence
                should exactly match with motif, _i.e._ Equal length and same
                positions for all `*`.

        """
        seq = ""
        struct = ""
        extra_bps = []
        pos = 0
        for i in range(len(motif)):
            c = motif[i]
            if c=="*":
                if sequence is not None and not sequence[i] == '*':
                    raise Exception("Motif and sequence are not compatible at position {}".format(i))
                extra_bps.append((pos + 1, pos + 2))
                seq += " & "
                struct += "(&)"
                pos += 2
            else:
                if sequence is not None:
                    w = sequence[i]
                else:
                    w = " "
                if w == '*':
                    raise Exception("Motif and sequence are not compatible at position {}".format(i))
                seq += w
                struct += c
                pos += 1
        seq = " " + seq + " "
        struct = "(" + struct + ")"
        self.sequence = seq
        self.structure = struct
        self.length = pos + 2
        extra_bps.append((0, self.length - 1))
        self.extra_bps = extra_bps

        self._init_features()
        # Default Bases Styles
        self.rootBasesStyle = BasesStyle(fill="#606060", outline="#FFFFFF",number="#FFFFFF")
        self.dummyBasesStyle = BasesStyle(fill="#DDDDDD", outline="#FFFFFF", number="#FFFFFF")

        self.default_color['baseNum'] = "#FFFFFF"
        self.params['bpStyle'] = 'simple'
        self.params['rotation'] = 180

    def _gen_input_cmd(self):
        return " -sequenceDBN \"{}\" -structureDBN \"{}\"".format(self.sequence, self.structure)

    def set_dummy_bases_style(self, style):
        """Set style for dummy bases. Argument is a [BasesStyle][varnaapi.BasesStyle] object.
        """
        if not isinstance(style, BasesStyle):
            raise Exception('The argument should be BasesStyle object')
        self.dummyBasesStyle = style

    def set_root_bases_style(self, style):
        """Set style for root bases. Argument is a [BasesStyle][varnaapi.BasesStyle] object.
        """
        if not isinstance(style, BasesStyle):
            raise Exception('The argument should be BasesStyle object')
        self.rootBasesStyle = style

    def savefig(self, output):
        dummybps = []
        for (i,j) in self.extra_bps:
            dummybps += [i, j]
            self.add_aux_BP(i=i, j=j, color="#DDDDDD", thickness=1)
        self.add_aux_BP(i=1, j=self.length-2, color="#000000", thickness=2)

        self.add_bases_style(self.rootBasesStyle, [1, self.length-2])
        self.add_bases_style(self.dummyBasesStyle, dummybps)
        super().savefig(output)
