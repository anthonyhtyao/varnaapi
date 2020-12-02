import re
import os
from string import ascii_lowercase, ascii_uppercase



VARNA_PATH="VARNAv3-93.jar"
HEX = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')
DEFAULT_COLOR_LIST = ['backbone', 'background', 'baseInner', 'baseName', 'baseNum',
                      'baseOutline', 'bp', 'gapsColor', 'nsBasesColor']

PARENTHESES_SYSTEMS = [
    ("(", ")"),
    ("[", "]"),
    ("<", ">"),
    ("{", "}")
] + [(c1, c2) for c1, c2 in zip(ascii_uppercase, ascii_lowercase)]
PARENTHESES_OPENING = [c1 for c1, c2 in PARENTHESES_SYSTEMS]
PARENTHESES_CLOSING = {c2: c1 for c1, c2 in PARENTHESES_SYSTEMS}

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
    def __init__(self, seq=None, structure=None):
        self.length = -1
        self.structure = []
        self.sequence = ""
        self.dbn = None
        self.aux_BPs = []
        self.highlight_regions = []
        self.algorithm = 'naview'
        self.default_color = {}
        self.title = None
        self.zoom = None

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
        if self.dbn is not None:
            self.dbn += "."*(self.length-len(self.structure))


    def add_aux_BP(self, i, j, edge5='wc', edge3='wc', stericity='cis', color='#0000FF', thickness=1.0):
        assert_valid_interval(self.length, i, j)
        if edge5 not in ['wc', 's', 'h']:
            raise Exception("edge5 should be one of wc, s, and h")
        if edge3 not in ['wc', 's', 'h']:
            raise Exception("edge3 should be one of wc, s, and h")
        if stericity not in ['cis', 'trans']:
            raise Exception("stericity should be either cis or trans")
        color = assert_hex_color(color)
        assert_is_number('thickness', thickness)

        self.aux_BPs.append((i, j, color, thickness, edge5, edge3, stericity))

    def add_highlight_region(self, i, j, fill="#9999FF", radius=15.0, outline="#3333FF"):
        assert_valid_interval(self.length, i, j)
        fill = assert_hex_color(fill)
        outline = assert_hex_color(outline)
        assert_is_number('radius', radius)

        self.highlight_regions.append((i, j, radius, fill, outline))

    def set_algorithm(self, algo):
        if algo not in ['line', 'circular', 'radiate', 'naview']:
            raise Exception("Sould be one of line, circular, radiate or naview")
        self.algorithm = algo


    def set_title(self, title, color='#808080', size=10):
        self.title = (title, color, size)

    def set_zoom_level(level):
        self.zoom = level

    def set_default_color(self, **kwargs):
        for key, value in kwargs.items():
            if not key in DEFAULT_COLOR_LIST:
                raise Exception("{} is not a valid keyword".format(key))
            assert_hex_color(value)
        self.default_color = kwargs

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

    def _gen_command(self):
        cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd".format(VARNA_PATH)
        cmd += " -sequenceDBN \"{}\" -structureDBN \"{}\" -o {}".format(self.sequence, self.format_structure(), self.output)
        cmd += " -algorithm {}".format(self.algorithm)

        # Command for defualt colors
        for key, color in self.default_color.items():
            if color is not None:
                cmd += " -{} {}".format(key, color)

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


        return cmd

    def savefig(self, output):
        self.output = output
        cmd = self._gen_command()
        print(cmd)
        os.popen(cmd).close()


    def __repr__(self):
        return repr((self.format_structure(),self.sequence))
