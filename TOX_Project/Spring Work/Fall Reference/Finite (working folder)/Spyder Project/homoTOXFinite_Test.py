import sys
import openmc
import math
import numpy as np
import pandas as pd
from pathlib import Path
import homoTOXMaterialModuleInput_Test as matMixMod
import matplotlib.pyplot as plt
import openmc.deplete


matMixOut = matMixMod.matMixFunInput()