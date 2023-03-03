from __future__ import annotations

import dataclasses
import json
from typing import List, Dict


@dataclasses.dataclass
class DACellGroup:
    id: str
    ncells: int
    lib_mean: float
    accessibilities: List[DAAccessibility] = None


@dataclasses.dataclass
class DAFeatureGroup:
    """
    """
    id: str
    proportion: float


@dataclasses.dataclass
class DAAccessibility:
    """
    """
    feature_group_id: str
    score: float

    def __post_init__(self):
        if 0 > self.score or self.score > 1:
            raise ValueError('accessibility score has to be between 0 and 1')


@dataclasses.dataclass
class DAConfig:
    cell_groups: Dict[str, DACellGroup]
    feature_groups: Dict[str, DAFeatureGroup]

    def __post_init__(self):
        proportions = [fg.proportion for fg in self.feature_groups.values()]
        for prop in proportions:
            if 1 < prop or prop < 0:
                raise ValueError('proportions cant be negative or grater than 1')
        if sum(proportions) > 1:
            raise ValueError('the total proportions need to sum up to maximum 1')

        # validate feature group ids
        for fg_id in [acc.feature_group_id for cg in self.cell_groups.values() for acc in cg.accessibilities]:
            if fg_id not in self.feature_groups:
                raise ValueError('feature group with id %s has not been defined' % fg_id)

        for acc in self.cell_groups.values():
            if len(acc.accessibilities) != len(self.feature_groups):
                raise ValueError('each cell group should have the same amount of accesssibilities '
                                 'as the number of feature groups')

    @classmethod
    def from_template(cls, file: str) -> Dict[str, DAConfig]:
        configs = {}
        with open(file, 'r') as f:
            data = json.loads(f.read())
            for simulation in data["simulations"]:
                # define cell groups
                cg_to_lib_mean = {lib_mean['id']: lib_mean['value'] for lib_mean in simulation['lib_means']}
                cgfg_to_score = {f'{score["cell_group"]}_{score["feature_group"]}': score['value']
                                 for score in simulation['scores']}
                cgs = {}
                for cg in data["cell_groups"]:
                    _id = cg['id']
                    cgs[_id] = DACellGroup(id=_id, ncells=cg['ncells'], lib_mean=cg_to_lib_mean[_id])

                # define feature groups
                fgs = {}
                for fg in data["feature_groups"]:
                    _id = fg['id']
                    fgs[_id] = DAFeatureGroup(id=_id, proportion=fg['proportion'])

                # define scores
                for cg_id, cg in cgs.items():
                    accs = []
                    for fg_id, fg in fgs.items():
                        score = cgfg_to_score[f'{cg_id}_{fg_id}']
                        accs.append(DAAccessibility(score=score, feature_group_id=fg_id))
                    cg.accessibilities = accs

                configs[simulation['id']] = DAConfig(cell_groups=cgs, feature_groups=fgs)
        return configs
