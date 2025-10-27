# -- import packages: ---------------------------------------------------------
import ABCParse
import pandas as pd

class GroupedData(ABCParse.ABCParse):
    def __init__(self, obs, groupby) -> None:

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
    def group1_idx(self) -> pd.Index:
        return self.group1.index

    @property
    def group2_idx(self) -> pd.Index:
        return self.group2.index

    def __call__(self, group1, group2) -> None:

        self.__update__(locals(), private=["group1", "group2"])
