import csv
from enum import Enum, auto
import math
import os
from typing import List, Tuple
import networkx as nx
from OCP.OCP import TopAbs
from OCP.OCP.BRepAdaptor import BRepAdaptor_Curve2d
from OCP.OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCP.OCP.BRepProj import BRepProj_Projection
from OCP.OCP.BRepTools import BRepTools_WireExplorer, BRepTools
from OCP.OCP.Geom import Geom_Plane
from OCP.OCP.TopExp import TopExp_Explorer
from OCP.OCP.TopoDS import TopoDS_Edge, TopoDS_Wire
from OCP.TopoDS import TopoDS
from ocp_vscode import show

from build123d import scale, Curve, SortBy, CenterOf
from build123d.build_enums import GeomType
from build123d.geometry import (
    TOLERANCE,
    Location,
    Plane,
    Pos,
    Rot,
    Vector,
    BoundBox,
    Axis,
)
from build123d.objects_curve import BSpline
from build123d.topology.composite import Compound
from build123d.topology.one_d import Edge, Wire, topo_explore_connected_faces
from build123d.topology.shape_core import ShapeList
from build123d.topology.three_d import Solid
from build123d.topology.two_d import Face, faces_are_tangent
from build123d.topology.zero_d import Vertex
from OCP.BRep import BRep_Tool


def unbend(
    bend_sequences: List[List[Face]],
    estimated_thickness: float,
    adj_graph: nx.Graph,
    seam_edges: ShapeList[Edge],
) -> ShapeList[Edge]:

    def flange_face_rotation(
        lcs: Plane,
        child_rot_seams: ShapeList[Edge],
        parent_face_bend_seams: ShapeList[Edge],
        bend_angle: float,
    ) -> Location:
        to_local = Location(lcs).inverse()
        to_world = Location(lcs)

        closest_child_rot_seam = child_rot_seams.sort_by_distance(lcs.origin)[0]
        ref_face_bend_seam = parent_face_bend_seams.sort_by_distance(lcs.origin)[0]

        child_seam_translation: Vector = (
            to_local * closest_child_rot_seam
        ).center() - (to_local * ref_face_bend_seam).center()
        transformation = Rot(bend_angle, 0, 0) * Pos(
            0, -child_seam_translation.Y, -child_seam_translation.Z
        )

        return to_world * transformation * to_local

    def bend_allowance_translation(
        lcs: Plane,
        bend_radius: float,
        bend_angle: float,
        thickness: float,
        k_factor: float | None = None,
    ) -> (Location, float):
        to_local = Location(lcs).inverse()
        to_world = Location(lcs)

        bac = BendAllowanceCalculator.read_file()
        if k_factor is None:
            k_factor = bac.get_k_factor(bend_radius, thickness)
        bend_allowance = (bend_radius + k_factor * thickness) * bend_angle

        bend_direction = 1 if bend_angle > 0 else -1

        return (
            to_world * Pos(0, -bend_direction * bend_allowance, 0) * to_local,
            bend_direction * bend_allowance,
        )

    flattened_edges = set()
    child: Face
    bend: Face
    parent: Face
    for unbend_path in bend_sequences:
        previous: Location = Location()
        for parent, bend, child in zip(
            unbend_path[0::2],  # Elements at indices 0, 2, 4, ...
            unbend_path[1::2],  # Elements at indices 1, 3, 5, ...
            unbend_path[2::2],  # Elements at indices 2, 4, 6, ...
        ):

            parent_bend_seam: Edge = adj_graph[parent][bend]["label"]
            child_bend_seam: Edge = adj_graph[bend][child]["label"]
            # Construct local coordinate system (border trihedron)
            local_origin = parent_bend_seam.end_point()
            local_y = parent.center() - parent_bend_seam.center()
            local_z = parent.normal_at(local_origin)
            local_x = local_y.cross(local_z)

            lcs = Plane(
                origin=local_origin,
                x_dir=local_x,
                z_dir=local_z,
            )

            signed_bend_angle: float = parent.normal_at().get_angle(child.normal_at())
            to_local = Location(lcs).inverse()

            rotation_direction = (
                (to_local * child).normal_at().cross((to_local * parent).normal_at())
            )

            if rotation_direction == Vector(0,0,0):
                # Sometimes we identify two faces with same normal sharing a bend, not necessary to unbend these
                flattened_edges.update(child.edges())
                continue


            bend_angle: float = (signed_bend_angle * rotation_direction).X
            current_flattened_edges = []
            rotation: Location = flange_face_rotation(
                lcs,
                ShapeList([child_bend_seam]),
                ShapeList([parent_bend_seam]),
                bend_angle,
            )
            bend_allowance_transform, bend_allowance_scaling = (
                bend_allowance_translation(
                    lcs, bend.radius, bend_angle * math.pi / 180, estimated_thickness
                )
            )

            unwrapped = previous * unwrap_cylinder(lcs, bend, bend_allowance_scaling)
            current_flattened_edges.extend(unwrapped)
            current: Location = bend_allowance_transform * rotation
            unbent_child_edges: ShapeList[Edge] = current * child.edges()
            current_flattened_edges.extend(previous * unbent_child_edges)
            current_flattened_edges.extend(previous * parent.edges())

            previous *= current
            flattened_edges.update(current_flattened_edges)
    return ShapeList(flattened_edges)

def uv_face(face: Face, projection_plane: Plane) -> Face:
    topods_outer_wire = BRepTools.OuterWire_s(face.wrapped)

    def uv_wire(wire: TopoDS_Wire) -> Wire:
        wire_builder = BRepBuilderAPI_MakeWire()
        edge_explorer = TopExp_Explorer(wire, TopAbs.TopAbs_EDGE)

        while edge_explorer.More():
            current_edge: TopoDS_Edge = TopoDS.Edge(edge_explorer.Current())
            curve_adaptor = BRepAdaptor_Curve2d(current_edge, face.wrapped)
            first, last = curve_adaptor.FirstParameter(), curve_adaptor.LastParameter()
            pcurve = BRep_Tool.CurveOnSurface_s(current_edge, face.wrapped, first, last)
            # TODO: check required types of projection plane wrapped - underlying surface?
            edge_builder = BRepBuilderAPI_MakeEdge(pcurve, Geom_Plane(projection_plane.wrapped), first, last)

            if edge_builder.IsDone():
                flatEdge: TopoDS_Edge = TopoDS.Edge(edge_builder.Edge())
                flatEdge.Orientation(current_edge.Orientation())
                wire_builder.Add(flatEdge)
            edge_explorer.Next()

        return Wire(wire_builder.Wire())

    outer_wire = uv_wire(topods_outer_wire)

    inner_wires = []
    wire_explorer = TopExp_Explorer(face.wrapped, TopAbs.TopAbs_WIRE)
    while wire_explorer.More():
        current_wire: TopoDS_Wire = TopoDS.Wire(wire_explorer.Current())
        if current_wire.IsSame(topods_outer_wire):
            wire_explorer.Next()
            continue

        inner_wires.append(uv_wire(current_wire))
        wire_explorer.Next()

    f = Face(outer_wire, inner_wires)
    return f

def unwrap_cylinder(lcs: Plane, bend: Face, bend_allowance: float) -> ShapeList[Edge]:
    to_local = Location(lcs).inverse()
    to_world = Location(lcs)


    u_min, u_max, v_min, v_max = bend._uv_bounds()

    # uw_face: Face = bend._uv_face(lcs, force_ordering=True)

    uw_lcs = Plane(
        origin=lcs.origin,
        x_dir=lcs.x_dir,
        z_dir=-lcs.z_dir if bend.is_circular_concave else lcs.z_dir,
    )

    #uw_face = uv_face(bend, uw_lcs)
    uw_face = bend._uv_face(uw_lcs, force_ordering=True)

    v_bend = to_local * Vertex(bend.position_at(0, 0)) - to_local * Vertex(bend.position_at(1, 0))
    v_uv = to_local * Vertex(uw_face.position_at(0, 0)) - to_local * Vertex(uw_face.position_at(1, 0))
    v_bend_xy = Vector(v_bend.X, v_bend.Y)
    v_uv_xy = Vector(v_uv.X, v_uv.Y)
    show(lcs, bend, uw_face)

    angle = v_bend_xy.get_signed_angle(v_uv_xy)

    show(lcs, bend, uw_face)

    bend_edges = bend.edges()
    uw_edges = uw_face.edges()

    show(bend_edges, uw_edges, lcs, bend_edges[3].start_point(), uw_edges[3].start_point())
    rot_corr = Rot(0,0,angle)

    uvs = [[0,0], [0,1], [1,0], [1,1]]
    uv_ref = sorted(
        map(
            lambda uv: (bend.position_at(uv[0], uv[1]), uv), uvs
        ),
        key=lambda pos_uv: (pos_uv[0] - lcs.origin).length
    )

    uw_face = to_world * rot_corr * to_local * uw_face

    show(bend, lcs, uw_face, Vertex(uv_ref[0][0]), Vertex(uw_face.position_at(uv_ref[0][1][0], uv_ref[0][1][1])))

    pos_corr = Pos(uv_ref[0][0] - uw_face.position_at(uv_ref[0][1][0], uv_ref[0][1][1]))

    correction = pos_corr

    uw_face: Face = correction * uw_face
    show(lcs, bend, uw_face)
    uw_edges = uw_face.edges()

    # Scale to account for projection factor and bend allowance
    allowance_factor: float = bend_allowance / (u_max - u_min)
    unwrapped_local = to_local * uw_edges
    scaled_edges = to_world * scale(
        unwrapped_local, (1, allowance_factor, 1), about=Vector(0, 0, 0)
    )
    show(lcs, bend, scaled_edges)
    return scaled_edges.edges()


def compute_bend_sequences(
    reference_face: Face, dfs_tree: nx.DiGraph, face_adjacency: nx.Graph
) -> Tuple[List[ShapeList[Face]], set[Edge]]:

    bend_sequences: List[ShapeList[Face]] = []
    seam_edges: set[Edge] = set()
    for u, v in dfs_tree.edges():
        seam_edges.add(face_adjacency[u][v]["label"])
        if dfs_tree.out_degree(v) == 0:
            # it's a leaf
            unfold_path = ShapeList(nx.shortest_path(dfs_tree, reference_face, v))
            is_valid_path = True
            for i in range(0, len(unfold_path) - 2, 2):
                first_flange: Face
                bend: Face
                second_flange: Face

                first_flange, bend, second_flange = (
                    unfold_path[i],
                    unfold_path[i + 1],
                    unfold_path[i + 2],
                )
                if not (
                    isinstance(first_flange.is_planar, Plane)
                    and (bend.is_circular_convex or bend.is_circular_concave)
                    and isinstance(second_flange.is_planar, Plane)
                ):
                    is_valid_path = False

            if is_valid_path:
                bend_sequences.append(unfold_path)
            else:
                raise RuntimeError(
                    f"Invalid pattern at indices {i}-{i + 2}: expected flange -> bend -> flange"
                )
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

    bend_sequences, seam_edges = compute_bend_sequences(
        reference_face, dfs_tree, tangent_faces_adjacacency_graph
    )
    estimated_thickness = estimate_thickness(solid_to_unfold, reference_face)
    flattened_edges = unbend(
        bend_sequences,
        estimated_thickness,
        tangent_faces_adjacacency_graph,
        ShapeList(seam_edges),
    )

    show(reference_face, flattened_edges, solid_to_unfold)

    return flattened_edges


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
                if len(connected_faces) == 2 and faces_are_tangent(
                    first=connected_faces[0],
                    second=connected_faces[1],
                    common_edge=f_edge,
                ):
                    face_1: Face = connected_faces[0]
                    face_2: Face = connected_faces[1]
                    if (
                        face_1.geom_type is GeomType.CYLINDER
                        and isinstance(face_2.is_planar, Plane)
                    ) or (
                        face_2.geom_type is GeomType.CYLINDER
                        and isinstance(face_1.is_planar, Plane)
                    ):
                        adjacent_faces_graph.add_node(
                            face_1, type=face_1.geom_type.__repr__()
                        )
                        adjacent_faces_graph.add_node(
                            face_2, type=face_2.geom_type.__repr__()
                        )
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

    def get_bend_allowance(
        self,
        radius: float,
        thickness: float,
        bend_angle: float,
    ) -> float:
        factor = self.get_k_factor(radius, thickness)
        bend_allowance = (radius + factor * thickness) * bend_angle
        return bend_allowance


def estimate_thickness(solid: Solid, ref_face: Face) -> float:
    # Get the normal of the reference face
    face_normal: Vector = ref_face.normal_at(ref_face.center())

    # Find all faces parallel to the reference face
    parallel_faces: ShapeList[Face] = ShapeList([
        f
        for f in solid.faces()
        if abs(f.normal_at(f.center()).dot(face_normal)) > (1 - TOLERANCE)
    ])

    # Filter for the opposite face (anti-parallel normal)
    opposite_faces = ShapeList([
        f
        for f in parallel_faces
        if f.normal_at(f.center()).dot(face_normal) < TOLERANCE
    ]).sort_by_distance(ref_face)

    if not opposite_faces:
        raise ValueError("No opposite face found for thickness estimation.")

    # Calculate the distance between the reference face and the opposite face
    thickness = ref_face.distance_to(opposite_faces[0])
    return thickness
