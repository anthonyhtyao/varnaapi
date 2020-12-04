import re
import os
from string import ascii_lowercase, ascii_uppercase



VARNA_PATH="VARNAv3-93.jar"
HEX = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')
BORDER = re.compile('^\d+x\d+$')
DEFAULT_COLOR_LIST = ['backbone', 'background', 'baseInner', 'baseName', 'baseNum',
                      'baseOutline', 'bp', 'gapsColor', 'nsBasesColor']
OPTIONS = ['autoHelices', 'autoInteriorLoops', 'autoTerminalLoops', 'drawBackbone', 'drawBases', 'drawNC', 'drawTertiary', 'fillBases', 'flat']

NUMERIC_OPTIONS = ['bpIncrement', 'periodNum', 'resolution', 'rotation', 'spaceBetweenBases', 'zoom']

PARENTHESES_SYSTEMS = [
    ("(", ")"),
    ("[", "]"),
    ("<", ">"),
    ("{", "}")
] + [(c1, c2) for c1, c2 in zip(ascii_uppercase, ascii_lowercase)]
PARENTHESES_OPENING = [c1 for c1, c2 in PARENTHESES_SYSTEMS]
PARENTHESES_CLOSING = {c2: c1 for c1, c2 in PARENTHESES_SYSTEMS}


class BasesStyle:
    def __init__(self, fill=None, outline=None, label=None, number=None):
        self.color = {}
        self.update(fill, outline, label, number)

    def update(self, fill=None, outline=None, label=None, number=None):
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

def bp_to_struct(bps):
    n = max([j for i, j in bps]) + 1
    ss = [-1 for i in range(n)]
    for i, j in bps:
        ss[i], ss[j] = j, i
    return ss

def parse_vienna(ss):
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
    """VARNA object
    """
    def __init__(self, seq=None, structure=None):
        """Test
        """
        self.length = -1
        self.structure = []
        self.sequence = ""
        self.dbn = None
        self.aux_BPs = []
        self.highlight_regions = []
        self.params = {'algorithm': "naview"}
        self.default_color = {}
        self.options = {}
        self.title = None
        self.bases_styles = {}

        if structure is not None:
            if isinstance(structure, list):
                if len(structure) > 0:
                    first = structure[0]
                    if len(first)==1:
                        self.structure = check_structure(structure)
                    elif len(first)==2:
                        self.structure = bp_to_struct(structure)
                    else:
                        raise Exception("Unrecognized structure format for %s"%(structure))
            # Dot-Bracket Notation
            elif isinstance(structure, str):
                self.structure = parse_vienna(structure)
                self.dbn = structure
            self.length = len(self.structure)
        if seq is not None:
            self.length = max(self.length,len(seq))
            self.sequence = seq
        # Now we know the length, let's extend the sequence and structure if necessary
        self.sequence += " "*(self.length-len(self.sequence))
        self.structure += [-1]*(self.length-len(self.structure))


    def add_aux_BP(self, i:int, j:int, edge5:str='wc', edge3:str='wc', stericity:str='cis', color:str='#0000FF', thickness:float=1.0):
        """Add an additional base pair `(i,j)`, possibly defining and using custom style

        Args:
            i: 5' position of base pair
            j: 3' position of base pair
            edge5: Edge 5' used for interaction in non-canonical base-pairs, as defined by the Leontis/Westhof classification of base-pairs. Admissible values are __wc__ (Watson/Crick edge), __h__ (Hoogsteen edge) and __s__ (Sugar edge).
            edge3: Edge 3' used for interaction in non-canonical base-pairs. Admissible values are __wc__, __h__ and __s__.
            stericity: Orientation of the strands. Admissible values are __cis__ and __trans__
            color: Base-pair color in Hex color codes
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

    def add_highlight_region(self, i:int, j:int, radius:float=15.0, fill:str="#9999FF", outline:str="#3333FF"):
        """Highlights a region by drawing a polygon of predefined radius,
        fill color and outline color around it.
        A region consists in an interval from base `i` to base `j`.

        Args:
            i: 5'-end of the highlight
            j: 3'-end of the highlight
            radius: Thickness of the highlight
            fill: The color used to fill the highlight
            outline: The color used to draw the line around the highlight
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
        if algo not in ['line', ' circular', 'radiate', 'naview']:
            raise Exception("Sould be one of line, circular, radiate or naview")
        self.params['algorithm'] = algo


    def set_title(self, title:str, color:str='#808080', size:int=10):
        """Set title displayed at the bottom of the panel with color and font size
        """
        self.title = (title, color, size)

    def set_zoom_level(self, level:float):
        """Defines the level of zoom and zoom increment used to display the RNA within this panel"""
        self.param['zoom'] = level

    def set_default_color(self, **kwargs):
        """Set default color used for different objects in the panel.

        Args:
            **kwargs: See below for the list of allowed keywords.
                      Value of a keyword is the default color, in Hex color codes, used for the related object.


        __Keywords:__

        | Name         | Object in panel                                         |
        |--------------|---------------------------------------------------------|
        | backbone     | phosphate-sugar backbone (aka skeleton) of the RNA      |
        | background   | background color used within the panel                  |
        | baseInner    | inner base color                                        |
        | baseName     | nucleotide name                                         |
        | baseNum      | base numbers                                            |
        | baseOutline  | outer base color                                        |
        | bp           | base-pair                                               |
        | nsBasesColor | non-standard bases (Anything but `A`, `C`, `G` or `U`)  |

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
            **kwargs: See below for detailed option lists.
                      Value of keyword is either `True` or `False`

        __Keywords:__

        Default value is set while [VARNA] object is created.

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
        for key, value in kwargs.items():
            if key not in OPTIONS:
                raise Exception("{} is not a valid keyword".format(key))
            if not isinstance(value, bool):
                raise Exception(key + "should be a boolean")
        self.options = kwargs

    def set_numeric_params(self, **kwargs):
        for key, value in kwargs.items():
            if key not in NUMERIC_OPTIONS:
                raise Exception("{} is not a valid keyword".format(key))
            assert_is_number(value)
            if key == "periodNum":
                self.params['periodNum'] = int(value)
            else:
                self.params[key] = float(value)

    def format_structure(self):
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
        """Set default style for base-pairs rendering, chosen among:

        | Label  | Description                                                                                            |
        |--------|--------------------------------------------------------------------------------------------------------|
        | none   | Base-pairs are not drawn, but can be implicitly seen from "ladders", _i.e_ helix structures            |
        | line   | A simple line is used to draw any base-pair, regardless of its type                                    |
        | rnaviz | A small square is drawn at equal distance of the two partners                                          |
        | lw     | Both canonical and non-canonical base-pairs are rendered according to the Leontis/Westhof nomenclature |

        __Note:__ `lw` is set by default
        """
        if style not in ['none', 'line', 'rnaviz', 'lw']:
            raise Exception('Should be one of none, line, rnaviz or lw')
        self.params['bpStyle'] = style

    def set_bp_increment(self, value):
        assert_is_number('value', value)
        self.params['bpIncrement'] = float(value)


    def add_bases_style(self, style, bases):
        if not isinstance(style, BasesStyle):
            raise Exception("style should be BasesStyle object")
        if len(bases) > 0:
            self.bases_styles[style] = self.bases_styles.get(style, set()).union({i+1 for i in bases})

    def _gen_command(self):
        """
        Return command to run VARNA
        """
        cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd".format(VARNA_PATH)
        cmd += " -sequenceDBN \"{}\" -structureDBN \"{}\" -o {}".format(self.sequence, self.format_structure(), self.output)

        # Command for defualt colors
        for key, color in self.default_color.items():
            if color is not None:
                cmd += " -{} {}".format(key, color)

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

        return cmd

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
