from pydantic import BaseModel
from typing import List


class Pathways(BaseModel):
    pathway_ids: List[str]
