import re
from typing import Any, Union

from lark.lexer import Token
from lark.visitors import Transformer


class BaseTransformer(Transformer):
    def NUM_ID(self, num: Token):
        return num.update(value=int(num))

    def TOKEN(self, tok: Token):
        return tok.update(value=tok.strip())

    def protein_id(self, children: list[Token]) -> Token:
        v = [x.value for x in children]
        return Token("protein_id", v)

    def ref_id(self, children: list[Token]) -> Token:
        v = [x.value for x in children]
        return Token("ref_id", v)

    @staticmethod
    def _fix_string(s: str) -> str:
        s = re.sub(r"([\[(]) ", r"\1", s)
        s = re.sub(r" ([\])])", r"\1", s)
        s = s.replace("\t", " ")
        return s


class GenericTreeTransformer(BaseTransformer):
    """Transform extracted values from bottom-up.

    Formats the tree into a dictionary with child token types as keys.
    """

    def description(self, children: list[Token]) -> Token:
        res = " ".join([x.value.strip() for x in children])
        res = self._fix_string(res)
        return Token("description", res)

    def content(self, children: list[Token]) -> dict[str, Union[str, list[int]]]:
        res = {x.type: x.value for x in children}
        return res

    def commentary(self, children: list[Token]) -> Token:
        return Token("commentary", children)

    def entry(self, children: list[Token]) -> dict[str, Any]:
        res = {x.type: x.value for x in children}
        return res


class ReactionTreeTransformer(GenericTreeTransformer):
    """Commentary in ``(...)`` are on substrates, and in ``|...|`` on products."""

    def reaction(self, children: list[Token]) -> Token:
        """Parse the reaction in the description node.

        There should be three parts:
            - The left-hand-side of the reaction, which is a list of chemical
                names, separated by ``<space>+<space>``.
            - The separator ``<space>=<space>``.
            - The right-hand-side of the reaction, which is a list of chemical
                names, separated by ``<space>+<space>``.
        """
        reaction = " ".join([x.value.strip() for x in children])
        reaction = self._fix_string(reaction)
        if " = " not in reaction:
            return Token("reaction_description", reaction)
        lhs, rhs = reaction.split(" = ")
        lhs = lhs.split(" + ")
        rhs = rhs.split(" + ")
        res = {"lhs": lhs, "rhs": rhs}
        return Token("reaction", res)

    def commentary(self, children: list[Token]) -> Token:
        return Token("commentary_substrate", children)

    def more_commentary(self, children: list[Token]) -> Token:
        return Token("commentary_product", children)

    def reversibility(self, children: list[Token]) -> Token:
        """``r`` for reversible, ``ir`` for irreversible, ``?`` for unknown."""
        x = children[0].value
        x = x.replace("{", "")
        x = x.replace("}", "")
        if not x:
            x = "?"
        return Token("reversibility", x)

    def entry(self, children: list[Token]) -> dict[str, Any]:
        res = {x.type: x.value for x in children}
        if "reversibility" in res:
            res["reaction"]["reversibility"] = res["reversibility"]
            del res["reversibility"]
        return res


class RefTreeTransformer(BaseTransformer):
    def ref_id(self, children: list[Token]) -> Token:
        return Token("ref_id", children[0].value)

    def citation(self, children: list[Token]) -> Token:
        res = " ".join([x.value.strip() for x in children])
        res = self._fix_string(res)
        return Token("citation", res)

    def pubmed(self, children: list[Token]) -> Token:
        res = children[0].value
        res = re.sub(r"\{Pubmed:(\d*)\}", r"\1", res)
        return Token("pubmed", res)

    def paper_stat(self, children: list[Token]) -> Token:
        res = children[0].split(",")
        return Token("paper_stat", res)

    def entry(self, children: list[Token]) -> dict[str, Any]:
        return {x.type: x.value for x in children}


class SpecificInfoTreeTransformer(GenericTreeTransformer):
    def substrate(self, children: list[Token]) -> Token:
        res = children[0].value.strip()
        return Token("substrate", res)


class CommentaryOnlyTreeTransformer(BaseTransformer):
    def description(self, children: list[Token]) -> Token:
        res = children[0]
        res = self._fix_string(res)
        res = res[1:-1]  # remove the brackets
        return Token("description", res)

    def entry(self, children: list[Token]) -> dict[str, Any]:
        res = {x.type: x.value for x in children}
        return res
