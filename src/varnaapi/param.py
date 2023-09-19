# Varna Parameter Setting

import yaml
from colour import Color

import varnaapi.settings

def _union_if_hack(*lst):
    """Return union set of two iterables if hack mode if enable, otherwise first one
    """
    if varnaapi.settings.CONFIG['hackmode']:
        return set.union(*[set(t) for t in lst])
    else:
        return lst[0]

# Define yaml-related function for Color object
def color_representer(dumper, data):
    return dumper.represent_scalar('!color', data.get_hex_l())


def color_constructor(loader, node):
    value = loader.construct_scalar(node)
    return Color(value)


yaml.add_representer(Color, color_representer)
yaml.add_constructor('!color', color_constructor)

# HEX = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')

#################
#               #
#     Color     #
#               #
#################

COLOR_LIST = ['backbone', 'background', 'baseInner', 'baseName', 'baseNum',
              'baseOutline', 'bp', 'gapsColor', 'nsBasesColor']
COLOR_LIST: list
"""Allowed options for basic color setting

    | Name         | Object in panel                                         | Default Color |
    |--------------|---------------------------------------------------------|---------|
    | backbone     | Phosphate-sugar backbone (aka skeleton) of the RNA      | <p style='color:#595959'>Cottonwood Gray #595959</p> |
    | background   | Background color used within the panel                  | <p style='background-color:gray;color:#FFFFFF'>White #FFFFFF</p> |
    | baseInner    | Inner base color                                        | <p style='background-color:gray;color:#F2F2F2'>White Smoke #F2F2F2</p> |
    | baseName     | Nucleotide name                                         | <p style='color:#000000'>Black #000000</p> |
    | baseNum      | Base numbers                                            | <p style='color:#3F3F3F'>Vibrant Black #3F3F3F</p> |
    | baseOutline  | Outer base color                                        | <p style='color:#595959'>Cottonwood Gray #595959</p> |
    | bp           | Base-pair color                                         | <p style='color:#0000FF'>Blue #0000FF</p> |
    | nsBasesColor | Non-standard bases (Anything but `A`, `C`, `G` or `U`)  | <p style='background-color:gray;color:#F2F2F2'>White Smoke #F2F2F2</p> |
"""

# TODO: Check gaps default color
COLOR_TYPE = {k: 'color' for k in COLOR_LIST}
COLOR_DEFAULT = {
    'backbone': Color ('#595959'),
    'background': Color('#FFFFFF'),
    'baseInner': Color('#F2F2F2'),
    'baseName': Color('#000000'),
    'baseNum': Color('#3F3F3F'),
    'baseOutline': Color('#595959'),
    'bp': Color('#0000FF'),
    'nsBasesColor': Color('#F2F2F2'),
    'gapsColor': Color('gray')
}


#################
#               #
#    Boolean    #
#               #
#################

BOOLEAN_OPTIONS = ['autoHelices', 'autoInteriorLoops', 'autoTerminalLoops', 'drawBackbone', 'drawBases', 'drawNC', 'drawTertiary', 'fillBases', 'flat']
BOOLEAN_OPTIONS: list
"""Boolean option list

    | Name              | Option                                                                                                                                        | Default |
    |-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|---------|
    | autoHelices       | Annotates each and every helix in the RNA with a unique `Hn` label                                                                            | False   |
    | autoInteriorLoops | Annotates each and every interior loop of the RNA with a unique `In` label                                                                    | False   |
    | autoTerminalLoops | Annotates each and every terminal loop of the RNA with a unique `Tn` label                                                                    | False   |
    | drawBackbone      | Backbone drawing                                                                                                                              | True    |
    | drawBases         | Displays the outline of a nucleotide base                                                                                                     | True    |
    | drawNC            | Displays non-canonical base-pairs                                                                                                             | True    |
    | drawTertiary      | Display of `non-planar` base-pairs, _i.e._ pseudoknots [^1]                                                                                   | True    |
    | fillBases         | Fill bases                                                                                                                                    | True    |
    | flat              | In `radiate` drawing mode, redraws the whole structure, aligning to a  baseline the base featured on the exterior loops (aka "dangling ends") | True    |

    [^1]: Since there is no canonical definition of pseudoknotted portions, a maximal planar subset is extracted from the input structure, defined to be the planar portion, and used as a scaffold for the drawing algorithms.

"""
BOOLEAN_TYPE = {k: bool for k in BOOLEAN_OPTIONS}
BOOLEAN_DEFAULT = {
    'autoHelices': False,
    'autoInteriorLoops': False,
    'autoTerminalLoops': False,
    'drawBackbone': True,
    'drawBases': True,
    'drawNC': True,
    'drawTertiary': True,
    'fillBases': True,
    'flat': True
}


#################
#               #
#    Numeric    #
#               #
#################

NUMERIC_TYPE = {'border': tuple, 'bpIncrement': float, 'periodNum': int, 'resolution': float, 'rotation': float, 'spaceBetweenBases': float, 'zoom': float}
NUMERIC_PARAMS = [k for k in NUMERIC_TYPE.keys()]
NUMERIC_PARAMS: list
"""Allowed numeric parameters

    | Name              | Type  | Description                                                                                                                                                                                                    | Default|
    |-------------------|-------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----|
    | border            | (int, int) | Sets the width and height of the panel border, <i>i.e.</i> the gap between the panel boundaries and those of the surface used to draw the RNA. Border setting is ignored if it's smaller than RNA draw. | N/A |
    | bpIncrement       | float | In linear drawing mode, defines the vertical increment used to separate two successive, nested base-pairs                                                                                                      | 0.65 |
    | periodNum         | int   | Sets the interval between two successive base numbers. More specifically, if `k` is the period, then the first and last bases  of the RNA are numbered, along with each base whose number is a multiple of `k` | 10 |
    | resolution        | float | Chooses the resolution of a bitmap PNG export, _i.e._ the multiplier in  the number of pixels in each dimension of the exported picture.                                                                           | 1 |
    | rotation          | float | Rotates the whole RNA of a certain angular increment                                                                                                                                                           | 0 |
    | spaceBetweenBases | float | Sets the distance between consecutive bases                                                                                                                                                                    | 1 |
    | zoom              | float | Defines the level of zoom and zoom increment used to display the RNA within this panel (min:0.5, max:60)                                                                                                                         | 1 |
"""

NUMERIC_DEFAULT = {
    'border': None,
    'bpIncrement': 0.65,
    'periodNum': 10,
    'resolution': 1,
    'rotation': 0,
    'spaceBetweenBases': 1,
    'zoom': 1
}


#################
#               #
#    Choices    #
#               #
#################

CHOICES_PARAMS = ['algorithm', 'bpStyle']
CHOICES_TYPE = {k: 'choices' for k in CHOICES_PARAMS}

BP_STYLES = ['none', 'simple', 'rnaviz', 'lw']
"""Allowed options for base-pair style (`bpStyle`), default value is `lw`

    | Label  | Description                                                                                                      |
    |--------|------------------------------------------------------------------------------------------------------------------|
    | none   | Base-pairs are not drawn, but can be implicitly seen from "ladders", _i.e_ helix structures                      |
    | simple | A simple line is used to draw any base-pair, regardless of its type                                              |
    | rnaviz | A small square is drawn at equal distance of the two partners                                                    |
    | lw     | Both canonical and non-canonical base-pairs are rendered according to the Leontis/Westhof nomenclature (Default) |


    Examples:
        >>> BasicDraw.update(bpStyle="simple")
"""

ALGORITHMS = ['line', 'circular', 'radiate', 'naview']
"""Allowed options for drawing algorithms (`algorithm`) are `line`, `circular`, `radiate`, and `naview`. The default value is `radiate`.

    Examples:
        >>> BasicDraw.update(algorithm="line")
"""

CHOICES_VALUE = {'algorithm': ALGORITHMS, 'bpStyle': BP_STYLES}
CHOICES_DEFAULT = {'algorithm': 'radiate', 'bpStyle': 'lw'}

# Union all parameters
PARAM_LIST = BOOLEAN_OPTIONS + COLOR_LIST + NUMERIC_PARAMS + CHOICES_PARAMS

PARAM_TYPE = BOOLEAN_TYPE.copy()
PARAM_TYPE.update(COLOR_TYPE)
PARAM_TYPE.update(NUMERIC_TYPE)
PARAM_TYPE.update(CHOICES_TYPE)

# Default value used in VARNA
_PARAM_DEFAULT = BOOLEAN_DEFAULT.copy()
_PARAM_DEFAULT.update(COLOR_DEFAULT)
_PARAM_DEFAULT.update(NUMERIC_DEFAULT)
_PARAM_DEFAULT.update(CHOICES_DEFAULT)

# Value from global loading
PARAM_DEFAULT = _PARAM_DEFAULT.copy()


def load_config(filename):
    """Load global style parameters

    Args:
        filename: path to file
    """
    global PARAM_DEFAULT
    with open(filename, 'r') as f:
        params = yaml.load(f, Loader=yaml.Loader)
    for par in _union_if_hack(PARAM_LIST, params.keys()):
        val = params.get(par, None)
        if val is not None:
            val = _params_type_check(par, val)
        PARAM_DEFAULT[par] = val


def _params_type_check(par, val):
    """Private parameter type check.
    In hack mode, check only color type
    """
    typ = PARAM_TYPE.get(par, None)
    if typ == 'color':
        try:
            val = Color(val)
        except ValueError as e:
            raise TypeError(str(e))
    # Turn off type chack in hack mode
    elif varnaapi.settings.CONFIG['hackmode']:
        return val
    elif typ is None:
        raise TypeError("Unknown parameter: {}".format(par))
    elif par == 'border':
        try:
            assert val is None or (isinstance(val[0], int) and isinstance(val[1], int))
        except AssertionError:
            raise TypeError("Value for border should be a pair of integers")
    elif typ == 'choices':
        if val not in CHOICES_VALUE[par]:
            raise TypeError('Value of {} should be one of {}'.format(par, CHOICES_VALUE[par]))
    else:
        try:
            val = typ(val)
        except:
            raise TypeError("The expected value type for {} is {} instead of {}".format(par, typ, type(val)))
    return val


class _DefaultObj:
    def __init__(self, **kwargs):
        self.params = list(kwargs.keys())
        self.default = kwargs.copy()
        self.values = {k: None for k in self.params}

    def _get_diff(self):
        res = {}
        for par, val in self.values.items():
            if not varnaapi.settings.CONFIG['hackmode']:
                assert par in self.params
            if val is not None and not val == self.default.get(par, None):
                if isinstance(val, Color):
                    res[par] = val.get_hex_l()
                else:
                    res[par] = val
        return res

    def _to_cmd(self):
        res = []
        for par, val in self._get_diff().items():
            res += ['-'+par, str(val)]
        return res


# Default Title
TITLE_DEFAULT = {'title': '', 'titleColor': Color('#000000'), 'titleSize': 19}

class _Title(_DefaultObj):
    def __init__(self, title, color, size, **kwargs):
        super().__init__(**TITLE_DEFAULT)
        try:
            assert not str(title) == ""
        except AssertionError:
            raise TypeError('Title cannot be empty string')
        self.values['title'] = str(title)
        self.values['titleColor'] = Color(color)
        self.values['titleSize'] = int(size)
        if varnaapi.settings.CONFIG['hackmode']:
            self.values.update(kwargs)


HIGHLIGHT_DEFAULT = {'radius': 16, 'fill': Color('#BCFFDD'), 'outline': Color('#6ED86E')}

class _Highlight(_DefaultObj):
    def __init__(self, radius, fill, outline, **kwargs):
        super().__init__(**HIGHLIGHT_DEFAULT)
        self.values['radius'] = float(radius)
        self.values['fill'] = Color(fill)
        self.values['outline'] = Color(outline)
        if varnaapi.settings.CONFIG['hackmode']:
            self.values.update(kwargs)

    def _to_cmd(self):
        return ','.join('{}={}'.format(k, v) for k, v in self._get_diff().items())


class BasesStyle(_DefaultObj):
    """Defines a custom base-style, to be applied later to a set of bases.
    A BasesStyle style contains colors used for different components of a base.
    BasesStyle is constructed from given colors for different components.

    Error:
        At least one argument should be given.

    Args:
        fill (color): color of inner part of base
        outline (color): color of outline of base
        label (color): base text (name) color
        number (color): base number color

    Examples:
        >>> style = BasesStyle(fill='#FF0000', outline='#00FF00')

    __See Also:__ [BasicDraw.add_bases_style][varnaapi.BasicDraw.add_bases_style]
    """
    def __init__(self, fill=None, outline=None, label=None, number=None, **kwargs):
        super().__init__(fill=COLOR_DEFAULT['baseInner'], outline=COLOR_DEFAULT['baseOutline'], label=COLOR_DEFAULT['baseName'], number=COLOR_DEFAULT['baseNum'])
        if varnaapi.settings.CONFIG['hackmode']:
            self._update(fill=fill, outline=outline, label=label, number=number, **kwargs)
        else:
            self._update(fill=fill, outline=outline, label=label, number=number)

    def _update(self, **kwargs):
        """Update component _colors.
        Same rule as [\_\_init\_\_][varnaapi.param.BasesStyle.__init__]
        """
        # if fill is None and outline is None and label is None and number is None:
        #     raise Exception("At least one should not be None")
        for par, val in kwargs.items():
            if val is not None:
                self.values[par] = Color(val)

    def _to_cmd(self, **kwargs):
        """Custom command generator for BasesStyle.
        Function takes the default bases color set by user in kwargs
        """
        for par, val in kwargs.items():
            self.default[par] = val
        return ",".join(k+"="+v for k, v in self._get_diff().items() if v is not None)


# For aux BP

BP_DEFAULT = {'edge5': 'wc', 'edge3': 'wc', 'stericity': 'cis', 'color': Color('blue'), 'thickness': 1}
BP_CHOICES = {'edge5': ['wc', 'h', 's'], 'edge3': ['wc', 'h', 's'], 'stericity': ['cis', 'trans']}

class _BPStyle(_DefaultObj):
    def __init__(self, **kwargs):
        super().__init__(**BP_DEFAULT)
        res = {}
        for key, val in kwargs.items():
            if key == 'color':
                res['color'] = Color(val)
            elif varnaapi.settings.CONFIG['hackmode']:
                res[key] = val
            elif key == 'thickness':
                res['thickness'] = float(val)
            elif key in BP_CHOICES:
                try:
                    assert val in BP_CHOICES[key]
                    res[key] = val
                except AssertionError:
                    raise TypeError('Value of {} should be one of {}'.format(key, BP_CHOICES[key]))
            else:
                raise TypeError('{} is not a valid keyword'.format(key))
        self.values = res


    def _to_cmd(self, color):
        self.default['color'] = color
        return ','.join('{}={}'.format(k, v) for k, v in self._get_diff().items())

#################
#               #
#  Annotation   #
#               #
#################

ANNOTATION_DEFAULT = {'type': 'L', 'color': Color('black'), 'size': 12}

class _Annotation(_DefaultObj):
    def __init__(self, text, aType, anchor, color="#000000", size=12, **kwargs):
        super().__init__(**ANNOTATION_DEFAULT)
        try:
            assert not str(text) == ""
        except AssertionError:
            raise TypeError('Text cannot be empty string')
        self.text = text
        self.values['type'] = aType
        self.values['color'] = Color(color)
        self.values['size'] = int(size)
        self.anchor = anchor

        if varnaapi.settings.CONFIG['hackmode']:
            self.values.update(kwargs)

    def _to_cmd(self):
        res = ["{}={}".format(k, v) for k, v in self._get_diff().items()]
        if isinstance(self.anchor, int):
            res.append("anchor={}".format(self.anchor))
        else:
            res += ["x={}".format(self.anchor[0]), "y={}".format(self.anchor[1])]
        return "{}:{}".format(self.text, ','.join(res))

class BaseAnnotation(_Annotation):
    """Annoation on a base.

    Args:
        text: Annotation caption
        anchor: Index of base to annotate
        color (color): Annotation color
        size (int): Font size
    """
    def __init__(self, text:str, anchor:int, color="#000000", size=12, **kwargs):
        super().__init__(text, 'B', int(anchor), color, size, **kwargs)

class LoopAnnotation(_Annotation):
    """Same as [BaseAnnotation][varnaapi.param.BaseAnnotation] but on a loop.
    Argument `anchor` can be index of any base in the loop of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=12, **kwargs):
        super().__init__(text, 'L', int(anchor), color, size,**kwargs)

class HelixAnnotation(_Annotation):
    """Same as [BaseAnnotation][varnaapi.param.BaseAnnotation] but on an helix.
    Argument `anchor` can be index of any base in the helix of interest.
    """
    def __init__(self, text, anchor, color="#000000", size=12, **kwargs):
        super().__init__(text, 'H', int(anchor), color, size, **kwargs)

class StaticAnnotation(_Annotation):
    """Annotation on a specified position in VARNA drawing.
    Unlike [BaseAnnotation][varnaapi.param.BaseAnnotation], argument `anchor` is omitted.
    However, arguments `x` and `y` are needed to specify annotation position.

    Danger:
        It is unrecommended to use static annotation unless you know what you're doing

    Args:
        x (int): x-coordinate of position
        y (int): y-ccordinate of position

    Examples:
        >>> sa = StaticAnnotation("Hello World", 100, 150, color="#FF0000")
    """
    def __init__(self, text, x, y, color="#000000", size=12, **kwargs):
        super().__init__(text, 'P', (int(x), int(y)),  color, size, **kwargs)

CHEM_DEFAULT = {'glyph': 'arrow', 'dir': 'in', 'intensity': 1, 'color': Color('#0000B2')}
CHEM_CHOICES = {'glyph': ['arrow', 'dot', 'pin', 'triangle'], 'dir': ['in', 'out']}

class _ChemProb(_DefaultObj):
    def __init__(self, **kwargs):
        super().__init__(**CHEM_DEFAULT)
        res = {}
        for key, val in kwargs.items():
            if key == 'color':
                res['color'] = Color(val)
            elif varnaapi.settings.CONFIG['hackmode']:
                res[key] = val
            elif key == 'intensity':
                res['intensity'] = float(val)
            elif key in CHEM_CHOICES:
                try:
                    assert val in CHEM_CHOICES[key]
                    res[key] = val
                except AssertionError:
                    raise TypeError('Value of {} should be one of {}'.format(key, CHEM_CHOICES[key]))
            else:
                raise TypeError('{} is not a valid keyword'.format(key))
        self.values = res


    def _to_cmd(self):
        return ','.join('{}={}'.format(k, v) for k, v in self._get_diff().items())

CM_DEFAULT = ["red", "blue", "green", "heat", "energy", "bw"]

class _ColorMap:
    def __init__(self, values, vMin, vMax, caption, style, **kwargs):
        self.values = list(map(float, values))
        if vMin is not None:
            vMin = float(vMin)
        self.vMin = vMin
        if vMax is not None:
            vMax = float(vMax)
        self.vMax = vMax
        self.caption = str(caption)
        if isinstance(style, str):
            if style not in CM_DEFAULT:
                raise ValueError('Value for colormap style should be one of {}'.format(CM_DEFAULT))
            self.style = style
        elif isinstance(style, dict):
            self.style = ';'.join("{}:{}".format(float(k), Color(v).get_hex_l()) for k,v in style.items())

        else:
            raise ValueError('Style should be either a string of {} or a dictionary'.format(CM_DEFAULT))
        if varnaapi.settings.CONFIG['hackmode']:
            self.kwargs = kwargs

    def _to_cmd(self):
        cmd = []
        cmd += ['-colorMap', ';'.join(map(str, self.values))]
        if self.caption != "":
            cmd += ['-colorMapCaption', self.caption]
        if self.vMin is not None:
            cmd += ['-colorMapMin', str(self.vMin)]
        if self.vMax is not None:
            cmd += ['-colorMapMax', str(self.vMax)]
        if self.style not in ['energy', '']:
            cmd += ['-colorMapStyle', self.style]
        if varnaapi.settings.CONFIG['hackmode']:
            for key, val in self.kwargs:
                cmd += ['-{}'.format(key), '{}'.format(val)]

        return cmd


#################
#               #
#     Main      #
#               #
#################


class _VarnaConfig:
    """Create default configuration for VARNA
    """
    def __init__(self):
        self._params = {key: None for key in PARAM_LIST}
        self._loaded_params = {key: None for key in PARAM_LIST}

    def _get_value(self, param):
        if self._params.get(param, None) is not None:
            return self._params.get(param, None)
        elif self._loaded_params.get(param, None) is not None:
            return self._loaded_params.get(param, None)
        elif PARAM_DEFAULT.get(param, None) is not None:
            return PARAM_DEFAULT.get(param, None)
        else:
            return _PARAM_DEFAULT.get(param, None)

    def _diff_param(self):
        """Get param name and value that is different than the default
        """
        params = {}
        for key in _union_if_hack(PARAM_LIST, self._params.keys(), self._loaded_params.keys()):
            val = self._get_value(key)
            default = _PARAM_DEFAULT.get(key, None)
            if val is not None and not val == default:
                params[key] = val
        return params

    def update(self, loaded=False, **kwargs):
        """Easy way to update params value.
        The list of keyward arguments can be found [here](config.md)
        """
        # Assert argument is in the parameter list and type check
        if loaded:
            params = self._loaded_params
        else:
            params = self._params
        for key, value in kwargs.items():
            if value is not None:
                try:
                    if not varnaapi.settings.CONFIG['hackmode']:
                        assert key in PARAM_LIST
                    value = _params_type_check(key, value)
                    params[key] = value
                except AssertionError:
                    print('{} is not a valid parameter name'.format(key))
                    print('A valid argument is one of', ', '.join(PARAM_LIST))
                except TypeError as e:
                    print(e)

    def get_params(self, complete:bool=False):
        """Get parameters with value in dictionary
        By default, only the parameters with value different than the default are returned.
        Set complete to True to get complete parameters.

        Args:
            complete: Return complete parameters.
        """
        if complete:
            param = self._params.copy()
        else:
            param = self._diff_param()
        return param

    # Individual update

    def set_algorithm(self, algo):
        """Set algorithm other than __radiate__ to draw secondary structure.
        Supported options are __line__, __circular__, __radiate__ and __naview__.
        """
        self.update(algorithm=algo)

    def set_bp_style(self, style):
        """Set default style for base-pairs rendering, chosen among [BP_STYLES][varnaapi.param.BP_STYLES]

        Note: `lw` is set by default

        Examples:
            >>> v = varnaapi.Structure()
            >>> v.set_bp_style("simple")
        """
        self.update(bpStyle=style)

    def _gen_param_cmd(self):
        """Return list of VARNA command generated from parameters
        """

        cmd = []
        for par, val in self.get_params().items():
            typ = PARAM_TYPE.get(par, None)
            if typ == 'color':
                cmd.append('-' + par)
                cmd.append(val.get_hex_l())
            elif par == 'border':
                cmd.append('-border')
                cmd.append('{}x{}'.format(*val))
            else:
                cmd.append('-' + par)
                cmd.append(str(val))
        return cmd

    def dump_param(self, filename):
        """Store style parameters into file

        Args:
            filename: path to the file
        """
        with open(filename, 'w') as f:
            yaml.dump(self.get_params(complete=True), f)

    def load_param(self, filename):
        """Load existing style parameters from file

        Args:
            filename: path to the file
        """
        with open(filename, 'r') as f:
            self.update(loaded=True, **yaml.load(f, Loader=yaml.Loader))


    # def set_zoom_level(self, level:float):
    #     """Defines the level of zoom and zoom increment used to display the RNA within this panel"""
    #     self.params['zoom'] = level

    # def set_default_color(self, **kwargs):
    #     """Set default color used for different objects in the panel.

    #     Args:
    #         **kwargs (dict): See [COLOR_LIST][varnaapi.COLOR_LIST] for the list of allowed keywords.
    #             Value of a keyword is the default color, in Hex color codes, used for the related object.

    #     Examples:
    #         >>> set_default_color({'backbone': '#000000', 'bp':'#FFFF00'})

    #     """
    #     for key, value in kwargs.items():
    #         if key not in COLOR_LIST:
    #             raise Exception("{} is not a valid keyword".format(key))
    #         assert_hex_color(value)
    #     self.default_color = kwargs

    # def toggle_options(self, **kwargs):
    #     """Enable or disable options of drawing

    #     Args:
    #         **kwargs (dict): See [OPTIONS][varnaapi.OPTIONS] for detailed option lists.
    #                   Value of keyword is either `True` or `False`


    #     """
    #     for key, value in kwargs.items():
    #         if key not in OPTIONS:
    #             raise Exception("{} is not a valid keyword".format(key))
    #         if not isinstance(value, bool):
    #             raise Exception(key + "should be a boolean")
    #     self.options = kwargs

    # def set_numeric_params(self, **kwargs):
    #     """Change value of numeric parameters in one function.
    #     This is equivalent to use setting function of each parameter,
    #     such as [set_bp_increment][varnaapi.VARNA.set_bp_increment].

    #     Args:
    #         **kwargs (dict): See [NUMERIC_PARAMS][varnaapi.NUMERIC_PARAMS] for allowed keywords.

    #     """
    #     for key, value in kwargs.items():
    #         if key not in NUMERIC_PARAMS:
    #             raise Exception("{} is not a valid keyword".format(key))
    #         assert_is_number(key, value)
    #         if key == "periodNum":
    #             self.params['periodNum'] = int(value)
    #         else:
    #             self.params[key] = float(value)

    # def set_border(self, border:str):
    #     """Sets the width and height of the panel border, _i.e._ the gap
    #     between the panel boundaries and those of the surface used to draw the RNA.
    #     Parameter `border` is in format `"wxh"` where `w` and `h` are width and height separated by `x`.

    #     Example:
    #         >>> set_border("20x30")
    #     """
    #     match = BORDER.search(border)
    #     if not match:
    #         raise Exception("border should be the format nxm where n and m are numbers")
    #     self.params['border'] = "\"{}\"".format(border)


    # def set_bp_increment(self, value:float):
    #     """In linear drawing mode, defines the vertical increment used to
    #     separate two successive, nested base-pairs.

    #     Example:
    #         >>> varna.set_bp_increment(1.2)
    #     """
    #     assert_is_number('value', value)
    #     self.params['bpIncrement'] = float(value)

