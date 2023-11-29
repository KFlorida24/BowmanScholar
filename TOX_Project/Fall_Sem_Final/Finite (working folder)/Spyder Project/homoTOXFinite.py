import sys
import openmc
import math
import numpy as np
import pandas as pd
from pathlib import Path
import homoTOXMaterialModuleInput as matMixMod
import matplotlib.pyplot as plt
import openmc.deplete
import xml.etree.ElementTree as et


matMixOut = matMixMod.matMixFunInput()