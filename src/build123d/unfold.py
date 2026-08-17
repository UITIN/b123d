import csv
from enum import Enum, auto
import math
import os
from typing import List, Tuple
import networkx as nx
from ocp_vscode import show

from build123d import scale
from build123d.build_enums import GeomType
from build123d.geometry import TOLERANCE, Location, Plane, Pos, Rot, Vector, BoundBox
from build123d.objects_curve import BSpline
from build123d.topology.composite import Compound
from build123d.topology.one_d import Edge, Wire, topo_explore_connected_faces
from build123d.topology.shape_core import ShapeList
from build123d.topology.three_d import Solid
from build123d.topology.two_d import Face, faces_are_tangent
from build123d.topology.zero_d import Vertex
from OCP.BRep import BRep_Tool


def unbend(bend_sequences: List[List[Face]], estimated_thickness: float, adj_graph: nx.Graph,
           seam_edges: ShapeList[Edge]) -> Compound:
    # todo,
    # return (t_r, r, t_c)
    # t_r - translation to match the seam edges of the flanges,
    # r - rotate flange by bend angle
    # t_c - bend compensation using k-factor (material dependent)

    # each of the

    def flange_face_rotation(lcs: Plane, child_rot_seams: ShapeList[Edge], parent_face_bend_seams: ShapeList[Edge],
                             bend_angle: float) -> Location:
        to_local = Location(lcs).inverse()
        to_world = Location(lcs)

        closest_child_rot_seam = child_rot_seams.sort_by_distance(lcs.origin)[0]
        ref_face_bend_seam = parent_face_bend_seams.sort_by_distance(lcs.origin)[0]

        child_seam_translation: Vector = (to_local * closest_child_rot_seam).center() - (
                    to_local * ref_face_bend_seam).center()
        transformation = Rot(0, bend_angle, 0) * Pos(-child_seam_translation.X, 0, -child_seam_translation.Z)

        return to_world * transformation * to_local

    def bend_allowance_translation(lcs: Plane, bend_radius: float, bend_angle: float, thickness: float,
                                   non_seam_cyl_edges: ShapeList[Edge] | None = None,
                                   k_factor: float | None = None) -> Location:
        to_local = Location(lcs).inverse()
        to_world = Location(lcs)

        bac = BendAllowanceCalculator.read_file()
        if k_factor is None:
            k_factor = bac.get_k_factor(bend_radius, thickness)
        bend_allowance = (bend_radius + k_factor * thickness) * bend_angle

        bend_direction = 1 if bend_angle > 0 else -1

        return to_world * Pos(-bend_direction * bend_allowance, 0, 0) * to_local

    child: Face
    bend: Face
    parent: Face
    for unbend_path in bend_sequences:
        previous: Location = Location()
        for parent, bend, child in zip(
                unbend_path[0::2],  # Elements at indices 0, 2, 4, ...
                unbend_path[1::2],  # Elements at indices 1, 3, 5, ...
                unbend_path[2::2]  # Elements at indices 2, 4, 6, ...
        ):
            # Your loop logic here
            # edges_to_transform.extend(filter(lambda e: e not in seam_edges, child.edges()))
            # edges_to_transform.extend(filter(lambda e: e not in seam_edges, bend.edges()))
            # not_seam.extend(filter(lambda e: e not in seam_edges, parent.edges()))

            parent_bend_seam: Edge = adj_graph[parent][bend]['label']
            child_bend_seam: Edge = adj_graph[bend][child]['label']

            # Construct local coordinate system (border trihedron)
            local_origin = parent_bend_seam.center()
            inside_vector = parent.center() - local_origin
            local_y = parent_bend_seam.tangent_at(0.5)
            local_z = parent.normal_at(local_origin)
            local_x = local_y.cross(local_z)

            if local_x.dot(inside_vector) < 0:
                local_x = -local_x


            lcs = Plane(
                origin=local_origin,
                x_dir=local_x,
                y_dir=local_y,
            )

            show(lcs, child, bend, parent)

            signed_bend_angle: float = parent.normal_at().get_signed_angle(
                child.normal_at()
            )
            to_local = Location(lcs).inverse()
            to_world = Location(lcs)

            rotation_direction = (to_local * child).normal_at().cross((to_local * parent).normal_at())
            bend_angle: float = (signed_bend_angle * rotation_direction).Y

            rotation: Location = flange_face_rotation(lcs, ShapeList([child_bend_seam]), ShapeList([parent_bend_seam]),
                                                      bend_angle)

            bend_allowance: Location = bend_allowance_translation(lcs, bend.radius, bend_angle * math.pi / 180,
                                                                  estimated_thickness)
            current: Location = previous * bend_allowance * rotation
            transformed_child: Face = current * child
            previous = current

            cyl_non_seam_edges = ShapeList(filter(lambda e: e not in seam_edges, bend.edges()))
            non_seam_child_edges = ShapeList(filter(lambda e: e not in seam_edges, transformed_child.edges()))

            non_seam_parent_edges = ShapeList(filter(lambda e: e not in seam_edges, parent.edges()))

            flattened_edges: ShapeList[Edge] = ShapeList(non_seam_child_edges + non_seam_parent_edges)
            cyl_bend_allowance_scaling: float = bend_allowance.position.Y / (bend.radius * bend_angle * (math.pi / 180))
            cyl_edge: Edge

            for cyl_edge in cyl_non_seam_edges:
                if cyl_edge.geom_type in [GeomType.LINE, GeomType.CIRCLE]:
                    start_vertex: Vector = transformed_child.vertices().sort_by_distance(cyl_edge)[0].position
                    end_vertex: Vector = parent.vertices().sort_by_distance(cyl_edge)[0].position
                    flattened_seam_edge: Edge = Edge.make_line(start_vertex, end_vertex)
                    flattened_edges.append(flattened_seam_edge)
                if cyl_edge.geom_type is GeomType.BSPLINE:
                    # spline_proj = project_bspline_to_plane(bend, cyl_edge, local_origin, local_x, local_z.cross(local_x), bend_allowance.position.Y, bend_angle)
                    # spline = to_world * to_local * flatten_edge_on_surface(lcs, bend, cyl_edge, bend_allowance.position.Y, bend_angle)
                    flattened_edges.append(flatten_bspline_2(lcs, bend, bend_allowance.position.Y, bend_angle))

            filtered = ShapeList(filter(lambda e: e not in seam_edges, flattened_edges))
            wire = Wire.combine(flattened_edges)
            x = 0


def flatten_bspline(lcs: Plane, bend: Face, bspline: Edge, bend_allowance: float, bend_angle: float):
    to_local = Location(lcs).inverse()
    to_world = Location(lcs)

    first, last = BRep_Tool.Range_s(bspline.wrapped, bend.wrapped)
    pcurve = BRep_Tool.Curve_s(bspline.wrapped, first, last)
    knots = list(map(lambda k: k, pcurve.KnotSequence()))
    poles = list(map(lambda p: Vector(p.X(), p.Y(), p.Z()), pcurve.Poles()))
    weights = pcurve.Weights() if pcurve.IsRational() else None

    is_periodic = pcurve.IsPeriodic()
    plane_normal: Vector = lcs.z_dir
    poles_flat_3d: List[Vector] = []
    pole: Vector
    for pole in poles:
        flat_pole_3d = pole - plane_normal * ((pole.dot(plane_normal)) / (plane_normal.dot(plane_normal)))
        poles_flat_3d.append(flat_pole_3d)
    flat_bspline_3d = BSpline(poles_flat_3d, knots, pcurve.Degree(), weights, pcurve.IsPeriodic())
    projection_factor: float = bspline.length / flat_bspline_3d.length
    allowance_factor: float = bend_allowance / (bend.radius * bend_angle * (math.pi / 180))

    return to_world * (Pos(0, projection_factor * allowance_factor, 0) * to_local * flat_bspline_3d)

def flatten_bspline_2(lcs: Plane, bend: Face, bend_allowance: float, bend_angle: float):
    to_local = Location(lcs).inverse()
    to_world = Location(lcs)
    show(bend, lcs)
    unwrapped_face = bend._uv_face(lcs)

    u_min, u_max, v_min, v_max = unwrapped_face._uv_bounds()

    show(unwrapped_face, bend, lcs)

    unwrapped_face = to_world * Pos(-u_max, 0, 0) * to_local * unwrapped_face
    show(unwrapped_face, bend, lcs)
    # Get the uv version of the face in parameter space: u from 0-2π
    u, v, w = tuple(BoundBox(unwrapped_face).size)
    # 6.283185307179586, 10.0

    # Scale u by radius
    unwrapped_face = scale(unwrapped_face, (bend.radius, 1, 1))
    x, y, z = tuple(BoundBox(unwrapped_face).size)
     # 31.415926735897937, 10.0000002
    show(unwrapped_face, bend, lcs)

    allowance_factor: float = bend_allowance / (bend.radius * bend_angle * (math.pi / 180))
    scaled = to_world * scale(unwrapped_face, (allowance_factor, 1, 1))


    return to_world * (Pos(0, allowance_factor, 0) * to_local)


def compute_bend_sequences(reference_face: Face, dfs_tree: nx.DiGraph, face_adjacency: nx.Graph) -> Tuple[
    List[ShapeList[Face]], set[Edge]]:

    bend_sequences: List[ShapeList[Face]] = []
    seam_edges: set[Edge] = set()
    for u, v in dfs_tree.edges():
        seam_edges.add(face_adjacency[u][v]['label'])
        if dfs_tree.out_degree(v) == 0:
            # it's a leaf
            unfold_path = ShapeList(nx.shortest_path(dfs_tree, reference_face, v))
            is_valid_path = True
            for i in range(0, len(unfold_path) - 2, 2):
                first_flange: Face
                bend: Face
                second_flange: Face

                first_flange, bend, second_flange = unfold_path[i], unfold_path[i + 1], unfold_path[i + 2]
                if not (isinstance(first_flange.is_planar, Plane) and (
                        bend.is_circular_convex or bend.is_circular_concave) and isinstance(second_flange.is_planar,
                                                                                            Plane)):
                    is_valid_path = False

            if is_valid_path:
                bend_sequences.append(unfold_path)
            else:
                raise RuntimeError(f"Invalid pattern at indices {i}-{i + 2}: expected flange -> bend -> flange")
    return bend_sequences, seam_edges


def unfold(solid_to_unfold: Solid, reference_face: Face, material: float) -> Solid:
    """Unfolds a solid given a reference face, on which plane we unfold the other faces

    Args:
        solid_to_unfold: Solid to unfold,
        reference_face: Face to which the unfold reference plane is constructed

    Returns:
        The unfolded part as a new Solid
    """

    tangent_faces_adjacacency_graph = build_graph(solid_to_unfold, reference_face)

    # depth first search tree with reference face as root
    dfs_tree = nx.dfs_tree(tangent_faces_adjacacency_graph, reference_face)

    (bend_sequences, seam_edges) = compute_bend_sequences(reference_face, dfs_tree, tangent_faces_adjacacency_graph)
    estimated_thickness = estimate_thickness(solid_to_unfold, reference_face)
    unbend(bend_sequences, estimated_thickness, tangent_faces_adjacacency_graph, ShapeList(seam_edges))

    unfolded_product = None
    bockar = None
    final_component = Compound([unfolded_product] + bockar)

    if final_component:
        show(final_component)

    return final_component


def build_graph(solid: Solid, root_face: Face) -> nx.Graph:
    adjacent_faces_graph = nx.Graph()
    for i, face in enumerate(solid.faces()):
        # if face.is_circular_concave or face.is_circular_convex:
        if face.geom_type is GeomType.CYLINDER:
            face_edges = face.edges()
            for f_edge in face_edges:
                connected_faces: ShapeList[Face] = ShapeList(
                    map(lambda f: Face(f), topo_explore_connected_faces(f_edge))
                )
                if len(connected_faces) == 2 and faces_are_tangent(first=connected_faces[0], second=connected_faces[1],
                                                                   common_edge=f_edge):
                    face_1: Face = connected_faces[0]
                    face_2: Face = connected_faces[1]
                    if (face_1.geom_type is GeomType.CYLINDER and isinstance(face_2.is_planar, Plane)) or \
                            (face_2.geom_type is GeomType.CYLINDER and isinstance(face_1.is_planar, Plane)):
                        adjacent_faces_graph.add_node(face_1, type=face_1.geom_type.__repr__())
                        adjacent_faces_graph.add_node(face_2, type=face_2.geom_type.__repr__())
                        adjacent_faces_graph.add_edge(
                            face_1,
                            face_2,
                            label=f_edge,
                        )
    # adjacent_faces_graph should have at least three connected subgraphs
    # (top side, bottom side, and sheet edge sides of the sheetmetal part).
    # We only care about the subgraph that includes the selected root face.

    for c in nx.connected_components(adjacent_faces_graph):
        if root_face in c:
            return adjacent_faces_graph.subgraph(c).copy()
    # If there is nothing tangent to the root face, return a graph with
    # one node and no edges.
    # This is useful for dxf/svg export of flat plates for manufacturing.
    single_face_graph = nx.Graph()
    single_face_graph.add_node(root_face)
    return single_face_graph


class BendAllowanceCalculator:
    class KFactorStandard(Enum):
        ANSI = auto()
        DIN = auto()

    def __init__(self) -> None:
        self.k_factor_standard = None
        self.radius_thickness_values = None
        self.k_factor = None

    @classmethod
    def read_file(cls):
        instance = cls()

        radius_thickness_list = []
        k_factor_list = []

        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, "k-factor.csv")

        with open(file_path, mode="r", encoding="utf-8") as file:
            reader = csv.reader(file)
            header = next(reader)
            a1 = header[0]
            b1 = header[1]
            r_t_header = "".join(c for c in a1 if c not in "' ").lower()
            if r_t_header != "radius/thickness":
                raise ValueError

            kf_header = "".join(c for c in b1 if c not in "' -()").lower()
            if kf_header == "kfactoransi":
                instance.k_factor_standard = cls.KFactorStandard.ANSI
            elif kf_header == "kfactordin":
                instance.k_factor_standard = cls.KFactorStandard.DIN
            else:
                raise ValueError

            for row in reader:
                if not row:
                    continue
                radius_thickness_list.append(float(row[0]))
                k_factor_list.append(float(row[1]))

        instance.radius_thickness_values = radius_thickness_list
        instance.k_factor_values = k_factor_list

        return instance

    def get_k_factor(self, radius, thickness):
        r_over_t = radius / thickness
        if r_over_t <= self.radius_thickness_values[0]:
            kf_val = self.k_factor_values[0]
        elif r_over_t >= self.radius_thickness_values[-1]:
            kf_val = self.k_factor_values[-1]
        else:
            i = 0
            while r_over_t <= self.radius_thickness_values[i]:
                i += 1
            kf1 = self.k_factor_values[i]
            kf2 = self.k_factor_values[i + 1]
            rt1 = self.radius_thickness_values[i]
            rt2 = self.radius_thickness_values[i + 1]
            kf_val = kf1 + (kf2 - kf1) * ((r_over_t - rt1) / (rt2 - rt1))
        return kf_val

    def get_bend_allowance(self, radius: float, thickness: float, bend_angle: float,
                           ) -> float:
        factor = self.get_k_factor(radius, thickness)
        bend_allowance = (radius + factor * thickness) * bend_angle
        return bend_allowance


def estimate_thickness(solid: Solid, ref_face: Face) -> float:
    # Get the normal of the reference face
    face_normal: Vector = ref_face.normal_at(
        ref_face.center()
    )

    # Find all faces parallel to the reference face
    parallel_faces: ShapeList[Face] = [
        f for f in solid.faces()
        if abs(f.normal_at(f.center()).dot(face_normal)) > (1 - TOLERANCE)
    ]

    # Filter for the opposite face (anti-parallel normal)
    opposite_faces = [
        f for f in parallel_faces
        if f.normal_at(f.center()).dot(face_normal) < TOLERANCE
    ]

    if not opposite_faces:
        raise ValueError("No opposite face found for thickness estimation.")

    # Calculate the distance between the reference face and the opposite face
    thickness = ref_face.distance_to(opposite_faces[0])
    return thickness