import unittest

from build123d import unfold
from build123d import import_step
from ocp_vscode import show


class UnfoldTests(unittest.TestCase):

    def test_build_graph(self):
        bend = import_step("/home/johan/Dev/b123d-fork/tests/hat.step")
        reF_face = bend.faces()[82]
        show(reF_face, bend, bend.faces())

        unfold(bend, reF_face, 1.0)

        # 41 på hat, 2 på U, 20 på folded cut, item1, item2, 3 på item3, 6 på item4, 5 på item5, 9 på item6
