# Varna Parameter Setting

from numbers import Real
from colour import Color


# HEX = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')

#################
#               #
#     Color     #
#               #
#################

COLOR_LIST = ['backbone', 'background', 'baseInner', 'baseName', 'baseNum',
              'baseOutline', 'bp', 'gapsColor', 'nsBasesColor']
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

# Default Bases Style

baseArgMap = {'fill': 'baseInner', 'outline': 'baseOutline', 'label': 'baseName', 'number': 'baseNumber'}


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
            label (Hex): base text (name) color
            number  (Hex): base number color

        Examples:
            >>> style = BasesStyle(fill='#FF0000', outline='#00FF00')
        """
        self._color = {}
        self.update(fill, outline, label, number)

    def update(self, fill=None, outline=None, label=None, number=None):
        """Update component _colors.
        Same rule as [\_\_init\_\_][varnaapi.BasesStyle.__init__]
        """
        # if fill is None and outline is None and label is None and number is None:
        #     raise Exception("At least one should not be None")
        if fill is not None:
            self._color["fill"] = Color(fill)
        if outline is not None:
            self._color["outline"] = Color(outline)
        if label is not None:
            self._color["label"] = Color(label)
        if number is not None:
            self._color["number"] = Color(number)

    def __str__(self):
        order = ['fill', 'outline', 'label', 'number']
        lst = ["{}={}".format(k, (self._color[k]).get_hex_l()) for k in order if k in self._color]
        return ",".join(lst)

    def cmd(self):
        order = ['fill', 'outline', 'label', 'number']
        lst = ["{}={}".format(k, (self._color[k]).get_hex_l()) for k in order if k in self._color and not self._color[k] == COLOR_DEFAULT[baseArgMap[k]]]
        return ",".join(lst)

#################
#               #
#    Boolean    #
#               #
#################

BOOLEAN_OPTIONS = ['autoHelices', 'autoInteriorLoops', 'autoTerminalLoops', 'drawBackbone', 'drawBases', 'drawNC', 'drawTertiary', 'fillBases', 'flat']
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

NUMERIC_TYPE = {'border': tuple, 'bpIncrement': int, 'periodNum': Real, 'resolution': Real, 'rotation': Real, 'spaceBetweenBases': Real, 'zoom': Real}
NUMERIC_PARAMS = [k for k in NUMERIC_TYPE.keys()]
"""Allowed numeric parameters
    | Label             | Type  | Description                                                                                                                                                                                                    | Default|
    |-------------------|-------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
    | border | (int, int) | Sets the width and height of the panel border, <i>i.e.</i> the gap between the panel boundaries and those of the surface used to draw the RNA. Border setting is ignored if it's smaller than RNA draw. | N/A |
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
"""Allowed options for base-pair style
    | Label  | Description                                                                                            |
    |--------|--------------------------------------------------------------------------------------------------------|
    | none   | Base-pairs are not drawn, but can be implicitly seen from "ladders", _i.e_ helix structures            |
    | simple | A simple line is used to draw any base-pair, regardless of its type                                    |
    | rnaviz | A small square is drawn at equal distance of the two partners                                          |
    | lw     | Both canonical and non-canonical base-pairs are rendered according to the Leontis/Westhof nomenclature |

    __See Also:__ [VARNA.set_bp_style][varnaapi.VARNA.set_bp_style]

"""

ALGORITHMS = ['line', 'circular', 'radiate', 'naview']

CHOICES_VALUE = {'algorithm': ALGORITHMS, 'bpStyle': BP_STYLES}
CHOICES_DEFAULT = {'algorithm': 'radiate', 'bpStyle': 'lw'}

# Union all parameters
PARAM_LIST = BOOLEAN_OPTIONS + COLOR_LIST + NUMERIC_PARAMS + CHOICES_PARAMS

PARAM_TYPE = BOOLEAN_TYPE.copy()
PARAM_TYPE.update(COLOR_TYPE)
PARAM_TYPE.update(NUMERIC_TYPE)
PARAM_TYPE.update(CHOICES_TYPE)

PARAM_DEFAULT = BOOLEAN_DEFAULT.copy()
PARAM_DEFAULT.update(COLOR_DEFAULT)
PARAM_DEFAULT.update(NUMERIC_DEFAULT)
PARAM_DEFAULT.update(CHOICES_DEFAULT)


def _params_type_check(**params):
    """Private parameter type check.
    """
    for par, val in params.items():
        typ = PARAM_TYPE[par]
        if typ == 'color':
            try:
                Color(val)
            except ValueError as e:
                raise TypeError(str(e))
        elif par == 'border':
            try:
                assert isinstance(val[0], int) and isinstance(type(val[1]), int)
            except AssertionError:
                raise TypeError("Value for border should be a pair of integers")
        elif typ == 'choices':
            if val not in CHOICES_VALUE[par]:
                raise TypeError('Value of {} should be one of {}'.format(par, CHOICES_VALUE[par]))
        else:
            if not isinstance(val, typ):
                raise TypeError("The expected value type for {} is {} instead of {}".format('par', typ, type(val)))




class VarnaConfig:
    """Create default configuration for VARNA
    """
    def __init__(self):
        self._params = PARAM_DEFAULT.copy()

    def _diff_param(self):
        """Get param name and value that is different than the default
        """
        params = {p: self._params[p] for p in PARAM_LIST if not self._params[p] == PARAM_DEFAULT[p]}
        if self._params['border'] is not None:
            params['border'] = "{}x{}".format(*self._params['border'])


        return params

    def update(self, **kwargs):
        """Easy way to update params value
        """
        # Assert argument is in the parameter list and type check
        for key, value in kwargs.items():
            try:
                assert key in PARAM_LIST
                _params_type_check(**{key: value})
                if PARAM_TYPE[key] == 'color':
                    self._params[key] = Color(value)
                else:
                    self._params[key] = value
            except AssertionError:
                print('{} is not a valid parameter name'.format(key))
                print('A valid argument is one of', ', '.join(PARAM_LIST))
            except TypeError as e:
                print(e)


    def get_params(self, complete=False):
        """Get parameters with value in dictionary
        By default, only the parameters with value different than the default are returned.
        Set complete to True to get complete parameters.

        Args:
            complete: Return complete parameters. Defaults to False
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
        """Set default style for base-pairs rendering, chosen among [BP_STYLES][varnaapi.BP_STYLES]

        __Note:__ `lw` is set by default

        Example:
            >>> varna.set_bp_style("simple")
        """
        self.update(bpStyle=style)

    def _gen_param_cmd(self):
        """Return list of VARNA command generated from parameters
        """

        cmd = []
        for par, val in self._diff_param().items():
            typ = PARAM_TYPE[par]
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

