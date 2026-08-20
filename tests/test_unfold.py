import unittest

from build123d import unfold
from build123d import import_step
import networkx as nx
from ocp_vscode import show


class UnfoldTests(unittest.TestCase):

    def test_build_graph(self):
        bend = import_step("/home/johan/Dev/b123d-fork/tests/item5.step")
        show(bend)
        reF_face = bend.faces()[5]

        unfold(bend, reF_face, 1.0)

        # 41 på hat, 2 på U, 20 på folded cut, item1, item2, 3 på item3, 6 på item4, 5 på item5, 9 på item6
