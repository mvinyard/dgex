

import pandas as pd


from . import utils

class GroupedData(utils.ABCParse):
    def __init__(self, obs, groupby):
        
        self.__parse__(locals())
        
    @property
    def _GROUPED(self):
        if not hasattr(self, "_grouped"):
            self._grouped = self.obs.groupby(self.groupby)
        return self._grouped
        
    @property
    def group1(self)-> pd.DataFrame:
        return self._GROUPED.get_group(self._group1)
        
    @property
    def group2(self)-> pd.DataFrame:
        return self._GROUPED.get_group(self._group2)
    
    @property
    def group1_idx(self):
        return self.group1.index
    
    @property
    def group2_idx(self):
        return self.group2.index
        
    def __call__(self, group1, group2):
        
        self.__update__(locals(), private=["group1", "group2"])