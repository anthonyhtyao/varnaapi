import re
import os
import abc
from string import ascii_lowercase, ascii_uppercase
from colour import Color
import subprocess

from param import VarnaConfig, BasesStyle, _Title, _Highlight, _Annotation

__version__ = '0.1.0'

_VARNA_PATH="VARNAv3-93.jar"


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


#################
#               #
#  Annotation   #
#               #
#################

# class _Annotation:
#     """Basic Annotation
#     """
#     def __init__(self, text, type, color="#000000", size=10):
#         self.text = text
#         self.type = type
#         self.color = color
#         self.size = size #: int: font size

#     def asdict(self):
#         return {'text': self.text, 'type': self.type, 'color': self.color,
#                 'size': self.size}

#     @abc.abstractmethod
#     def to_cmd(self):
#         pass


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


class BaseAnnotation(_Annotation):
    def __init__(self, text:str, anchor:int, color="#000000", size=12):
        """Annoation on a base.

        Args:
            text: Annotation caption
            anchor: Index of base to annotate
            color (Hex): Annotation color
            size (int): Font size
        """
        super().__init__(text, 'B', int(anchor), color, size)


class LoopAnnotation(_Annotation):
    """Same as [BaseAnnotation][varnaapi.BaseAnnotation] but on a loop.
    Argument `anchor` can be index of any base in the loop of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=12):
        super().__init__(text, 'L', int(anchor), color, size)


class HelixAnnotation(_ObjectAnnotation):
    """Same as [BaseAnnotation][varnaapi.BaseAnnotation] but on an helix.
    Argument `anchor` can be index of any base in the helix of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=12):
        super().__init__(text, 'H', int(anchor), color, size)


class StaticAnnotation(_Annotation):
    def __init__(self, text, x, y, color="#000000", size=12):
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
        super().__init__(text, 'P', (int(x), int(y)),  color, size)

# def is_hex_color(color):
#     match = HEX.search(color)
#     if match:
#         return True
#     return False

# def assert_hex_color(color):
#     if not is_hex_color(color):
#         raise Exception("{} is not a valid Hex color code".format(color))
#     return color.upper()


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

class VARNA(VarnaConfig):
    def __init__(self, **kwargs):
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
        super().__init__()
        self._init_features()
        self._read_input(**kwargs)

    def _read_input(self, seq=None, structure=None):
        self.length = -1
        self.structure = []
        self.sequence = ""

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
        self._title = None
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

    def add_highlight_region(self, i:int, j:int, radius:float=16, fill="#BCFFDD", outline="#6ED86E"):
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

        self.highlight_regions.append((i+1, j+1, _Highlight(radius, fill, outline)))

    def set_title(self, title:str, color='#000000', size:int=19):
        """Set title displayed at the bottom of the panel with color and font size
        """
        self._title = _Title(title, color, size)


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
        cmd = ['java', '-cp', _VARNA_PATH, 'fr.orsay.lri.varna.applications.VARNAcmd']
        # cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd".format(_VARNA_PATH)

        cmd += self._gen_input_cmd()

        cmd += ['-o', self.output]
        # cmd += " -o {}".format(self.output)

        cmd += self._gen_param_cmd()

        # Title cmd
        if self._title is not None:
            cmd += self._title.to_cmd()

        # Aux Base pairs
        # if len(self.aux_BPs) > 0:
        #     auxbps = ["({},{}):color={},thickness={},edge5={},edge3={},stericity={}"
        #               .format(*t) for t in self.aux_BPs]
        #     cmd += " -auxBPs \"{}\"".format(";".join(auxbps))

        # Highlight Region
        if len(self.highlight_regions) > 0:
            res = []
            for item in self.highlight_regions:
                s = "{}-{}".format(item[0], item[1])
                setting = item[2].to_cmd()
                if not setting == "":
                    s += ":"+setting
                res.append(s)
            cmd += ['-highlightRegion', ';'.join(res)]

        # BasesStyles
        styles = {'fill': 'baseInner', 'outline': 'baseOutline', 'label': 'baseName', 'number': 'baseNum'}
        styles_dafault = {v: self.get_params().get(v) for v in styles.values() if v in self.get_params()}
        for ind, (style, bases) in enumerate(self.bases_styles.items()):
            s = style.to_cmd(**styles_dafault)
            if not s == "":
                cmd += ["-basesStyle{}".format(ind+1), s]
                cmd += ["-applyBasesStyle{}on".format(ind+1), ','.join(map(str, bases))]

        # Annotations
        if len(self.annotations) > 0:
            cmd += ["-annotations", ';'.join([t.to_cmd() for t in self.annotations])]

        return cmd

    def _gen_input_cmd(self):
        return ['-sequenceDBN', self.sequence, '-structureDBN', self.format_structure()]
        # return " -sequenceDBN \"{}\" -structureDBN \"{}\"".format(self.sequence, self.format_structure())

    def savefig(self, output):
        """
        Call VARNA to draw and store the paint in output
        """
        self.output = output
        cmd = self._gen_command()
        print(cmd)
        # os.popen(cmd).close()
        subprocess.run(cmd)

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
