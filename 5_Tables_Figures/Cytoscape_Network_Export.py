#!/usr/bin/env python3
"""Build auditable Cytoscape inputs from the saved Model 1 pathway analysis.

This script does not rerun enrichment or clustering. Node attributes come from
the silhouette-selected pathway-profile clusters, and edges come from the
prespecified complete-membership overlap graph (overlap coefficient >= 0.50).
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "tables"
CYTOSCAPE = PROJECT / "cytoscape"

MEMBERSHIP = TABLES / "pathway_profile_cluster_membership_model1.csv"
REPRESENTATIVES = TABLES / "pathway_profile_cluster_representatives_model1.csv"
DRIVERS = TABLES / "pathway_profile_cluster_drivers_model1.csv"
EDGES = TABLES / "igraph_full_membership_edges_model1.csv"

NODE_OUTPUT = CYTOSCAPE / "pathway_nodes_model1.csv"
EDGE_OUTPUT = CYTOSCAPE / "pathway_edges_model1.csv"
FULL_CYJS_OUTPUT = CYTOSCAPE / "pathway_network_full_model1.cyjs"
DISPLAY_CYJS_OUTPUT = CYTOSCAPE / "pathway_network_display_model1.cyjs"
MANIFEST_OUTPUT = CYTOSCAPE / "pathway_network_manifest_model1.json"

REPRESENTATIVE_LABELS = {
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION": "Asparagine N-linked glycosylation",
    "REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS": "Growth factor receptor signaling",
    "REACTOME_CA2_PATHWAY": "Ca2+ pathway",
    "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_BROWN_AND_BEIGE_ADIPOCYTE_DIFFERENTIATION": "Brown/beige adipocyte differentiation",
    "HALLMARK_MYC_TARGETS_V2": "MYC targets V2",
    "REACTOME_TRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL": "Nuclear/cytosolic tRNA modification",
    "REACTOME_HCMV_EARLY_EVENTS": "HCMV early events",
    "REACTOME_SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION": "SIRT1 regulation of rRNA expression",
    "REACTOME_INTERFERON_SIGNALING": "Interferon signaling",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING": "Interferon-alpha/beta signaling",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    return float(value) if value not in ("", None) else math.nan


def as_int(value: str) -> int:
    return int(float(value)) if value not in ("", None) else 0


def clean_pathway_label(pathway: str) -> str:
    label = re.sub(r"^(HALLMARK|REACTOME)_", "", pathway)
    label = label.replace("_", " ").lower().title()
    replacements = {
        "Dna": "DNA",
        "Rna": "RNA",
        "Rrna": "rRNA",
        "Mrna": "mRNA",
        "Trna": "tRNA",
        "Nfkb": "NF-kB",
        "Tnfa": "TNF-alpha",
        "Myc": "MYC",
        "Mtorc1": "mTORC1",
        "Jak": "JAK",
        "Stat": "STAT",
        "Il ": "IL ",
        "Interferon Alpha Beta": "Interferon-alpha/beta",
        "Ca2 ": "Ca2+ ",
        "Hcmv": "HCMV",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
    return label


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def cytoscape_payload(nodes: list[dict[str, object]], edges: list[dict[str, object]], name: str) -> dict[str, object]:
    return {
        "data": {
            "name": name,
            "description": (
                "Saved Model 1 FGSEA pathway-profile clusters. Node color denotes the "
                "three-contrast profile cluster; edges denote gene-set overlap."
            ),
        },
        "elements": {
            "nodes": [{"data": node} for node in nodes],
            "edges": [{"data": edge} for edge in edges],
        },
    }


def main() -> None:
    CYTOSCAPE.mkdir(parents=True, exist_ok=True)
    membership = read_csv(MEMBERSHIP)
    representatives = read_csv(REPRESENTATIVES)
    drivers = read_csv(DRIVERS)
    source_edges = read_csv(EDGES)

    graph_degree: dict[str, int] = defaultdict(int)
    graph_strength: dict[str, float] = defaultdict(float)
    for edge in source_edges:
        weight = as_float(edge["weight"])
        for endpoint in (edge["from"], edge["to"]):
            graph_degree[endpoint] += 1
            graph_strength[endpoint] += weight

    representative_meta = {
        row["pathway"]: {
            "representative_role": row["representative_role"],
            "representative_rank": as_int(row["representative_rank"]),
            "selection_basis": row["selection_basis"],
        }
        for row in representatives
    }
    driver_by_cluster: dict[str, list[str]] = defaultdict(list)
    for row in sorted(drivers, key=lambda x: (as_int(x["cluster_number"]), as_int(x["driver_rank"]))):
        driver_by_cluster[row["cluster_id"]].append(row["gene"])

    nodes: list[dict[str, object]] = []
    for row in membership:
        pathway = row["pathway"]
        rep = representative_meta.get(pathway, {})
        role = str(rep.get("representative_role", "Not selected"))
        short_label = clean_pathway_label(pathway)
        nodes.append(
            {
                "id": pathway,
                "name": pathway,
                "pathway_label": short_label,
                "plot_label": REPRESENTATIVE_LABELS.get(pathway, short_label) if role != "Not selected" else "",
                "pathway_source": row["source"],
                "cluster_id": row["cluster_id"],
                "cluster_number": as_int(row["cluster_number"]),
                "newly_diagnosed_NES": as_float(row["DisCOPDAI_1vs0"]),
                "potentially_excluded_NES": as_float(row["DisCOPDAI_2vs0"]),
                "established_COPD_NES": as_float(row["DisCOPDAI_3vs0"]),
                "significant_contrasts": as_int(row["significant_contrasts"]),
                "minimum_FDR": as_float(row["minimum_fdr"]),
                "negative_log10_minimum_FDR": min(-math.log10(max(as_float(row["minimum_fdr"]), 1e-300)), 20),
                "aggregate_evidence": as_float(row["aggregate_evidence"]),
                "centroid_distance": as_float(row["profile_distance_to_centroid"]),
                "overlap_graph_degree": graph_degree[pathway],
                "overlap_graph_strength": graph_strength[pathway],
                "representative_role": role,
                "representative_rank": rep.get("representative_rank", 0),
                "selection_basis": rep.get("selection_basis", ""),
                "leading_edge_drivers": "; ".join(driver_by_cluster[row["cluster_id"]]),
            }
        )

    node_ids = {str(node["id"]) for node in nodes}
    edges: list[dict[str, object]] = []
    for index, row in enumerate(source_edges, start=1):
        source = row["from"]
        target = row["to"]
        if source not in node_ids or target not in node_ids:
            raise ValueError(f"Edge endpoint missing from membership table: {source} -- {target}")
        weight = as_float(row["weight"])
        jaccard = as_float(row["jaccard"])
        edge = {
            "id": f"e{index:04d}",
            "source": source,
            "target": target,
            "interaction": "gene-set overlap",
            "overlap_coefficient": weight,
            "jaccard_index": jaccard,
            "shared_gene_count": as_int(row["intersection_n"]),
            "display_edge": False,
        }
        edges.append(edge)

    # The manuscript view is balanced across the five profile clusters. It
    # always includes both representatives, then adds the most connected
    # pathways until each cluster has four nodes. This is a display rule, not
    # a new inferential pathway-selection procedure.
    adjacency: dict[str, dict[str, float]] = defaultdict(dict)
    for edge in edges:
        source, target = str(edge["source"]), str(edge["target"])
        weight = float(edge["overlap_coefficient"])
        adjacency[source][target] = weight
        adjacency[target][source] = weight

    display_node_ids: set[str] = set(representative_meta)
    for _ in range(2):
        for cluster_number in range(1, 6):
            candidates = [
                node for node in nodes
                if node["cluster_number"] == cluster_number and str(node["id"]) not in display_node_ids
            ]
            candidates.sort(
                key=lambda node: (
                    -sum(neighbor in display_node_ids for neighbor in adjacency[str(node["id"])]),
                    -sum(
                        weight for neighbor, weight in adjacency[str(node["id"])].items()
                        if neighbor in display_node_ids
                    ),
                    -int(node["overlap_graph_degree"]),
                    -float(node["overlap_graph_strength"]),
                    -int(node["significant_contrasts"]),
                    -float(node["aggregate_evidence"]),
                    str(node["name"]),
                )
            )
            display_node_ids.add(str(candidates[0]["id"]))

    induced_degree = {
        pathway: sum(neighbor in display_node_ids for neighbor in adjacency[pathway])
        for pathway in display_node_ids
    }
    display_node_ids = {
        pathway for pathway in display_node_ids
        if pathway in representative_meta or induced_degree[pathway] > 0
    }

    for node in nodes:
        selected = str(node["id"]) in display_node_ids
        node["display_selected"] = selected
        node["display_selection_reason"] = (
            str(node["representative_role"])
            if str(node["representative_role"]) != "Not selected"
            else ("Cluster overlap context" if selected else "Full network only")
        )

    display_edges = []
    for edge in edges:
        selected = str(edge["source"]) in display_node_ids and str(edge["target"]) in display_node_ids
        edge["display_edge"] = selected
        if selected:
            display_edges.append(edge)
    display_nodes = [node for node in nodes if bool(node["display_selected"])]

    node_fields = list(nodes[0].keys())
    edge_fields = list(edges[0].keys())
    write_rows(NODE_OUTPUT, nodes, node_fields)
    write_rows(EDGE_OUTPUT, edges, edge_fields)

    full_payload = cytoscape_payload(nodes, edges, "COPD Schema pathway network - full overlap graph")
    display_payload = cytoscape_payload(display_nodes, display_edges, "COPD Schema pathway network - manuscript view")
    FULL_CYJS_OUTPUT.write_text(json.dumps(full_payload, indent=2, allow_nan=False), encoding="utf-8")
    DISPLAY_CYJS_OUTPUT.write_text(json.dumps(display_payload, indent=2, allow_nan=False), encoding="utf-8")

    manifest = {
        "analysis": "Saved Model 1 pathway results; no enrichment or clustering recomputed",
        "node_source": str(MEMBERSHIP.relative_to(PROJECT)),
        "edge_source": str(EDGES.relative_to(PROJECT)),
        "nodes": len(nodes),
        "full_edges": len(edges),
        "display_nodes": len(display_nodes),
        "display_edges": len(display_edges),
        "edge_definition": "Overlap coefficient >= 0.50 in the saved complete-membership graph",
        "display_node_rule": (
            "Both representatives plus connected overlap context, with up to "
            "4 pathways per profile cluster; unlabeled nonrepresentative isolates omitted"
        ),
        "display_edge_rule": "All saved overlap edges induced by the display pathways",
        "label_rule": "Labels shown only for centroid-nearest and central-evidence representatives",
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

