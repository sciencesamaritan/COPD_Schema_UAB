#!/usr/bin/env python3
"""Render the saved pathway network in Cytoscape through CyREST.

Cytoscape Desktop must be running with CyREST listening on localhost:1234.
The script creates a fresh session, imports the display CYJS, applies a visual
style, uses a weighted force-directed layout, exports vector/raster figures,
and saves the editable Cytoscape session.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "cytoscape" / "pathway_network_display_model1.cyjs"
FIGURES = PROJECT / "figures"
SESSION = PROJECT / "cytoscape" / "pathway_network_model1.cys"
BASE = "http://localhost:1234/v1"
STYLE_NAME = "COPD Schema pathway clusters"


def request(path: str, *, method: str = "GET", payload: object | None = None) -> object:
    body = None if payload is None else json.dumps(payload).encode("utf-8")
    req = Request(
        f"{BASE}{path}",
        data=body,
        method=method,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    try:
        with urlopen(req, timeout=180) as response:
            raw = response.read().decode("utf-8")
    except (HTTPError, URLError) as error:
        detail = getattr(error, "read", lambda: b"")().decode("utf-8", errors="replace")
        raise RuntimeError(f"CyREST {method} {path} failed: {error}\n{detail}") from error
    return json.loads(raw) if raw.strip().startswith(("{", "[")) else raw


def command(namespace: str, name: str, **parameters: object) -> object:
    query = urlencode({key: str(value) for key, value in parameters.items()})
    suffix = f"?{query}" if query else ""
    return request(f"/commands/{quote(namespace)}/{quote(name)}{suffix}", method="POST")


def download(path: str, destination: Path) -> None:
    req = Request(f"{BASE}{path}", method="GET")
    try:
        with urlopen(req, timeout=180) as response:
            destination.write_bytes(response.read())
    except (HTTPError, URLError) as error:
        detail = getattr(error, "read", lambda: b"")().decode("utf-8", errors="replace")
        raise RuntimeError(f"CyREST GET {path} failed: {error}\n{detail}") from error


def style_definition() -> dict[str, object]:
    cluster_colors = {
        "C1": "#4477AA",
        "C2": "#CCBB44",
        "C3": "#AA3377",
        "C4": "#EE7733",
        "C5": "#228833",
    }
    return {
        "title": STYLE_NAME,
        "defaults": [
            {"visualProperty": "NETWORK_BACKGROUND_PAINT", "value": "#FFFFFF"},
            {"visualProperty": "NODE_FILL_COLOR", "value": "#B8BEC6"},
            {"visualProperty": "NODE_BORDER_PAINT", "value": "#222222"},
            {"visualProperty": "NODE_BORDER_WIDTH", "value": "0.8"},
            {"visualProperty": "NODE_SIZE", "value": "28.0"},
            {"visualProperty": "NODE_SHAPE", "value": "ELLIPSE"},
            {"visualProperty": "NODE_LABEL_COLOR", "value": "#111111"},
            {"visualProperty": "NODE_LABEL_FONT_FACE", "value": "Arial,plain,10"},
            {"visualProperty": "NODE_LABEL_WIDTH", "value": "160.0"},
            {"visualProperty": "NODE_LABEL_POSITION", "value": "S,N,c,0.00,4.00"},
            {"visualProperty": "EDGE_STROKE_UNSELECTED_PAINT", "value": "#A7ADB4"},
            {"visualProperty": "EDGE_TRANSPARENCY", "value": "75"},
            {"visualProperty": "EDGE_WIDTH", "value": "0.8"},
            {"visualProperty": "EDGE_LINE_TYPE", "value": "SOLID"},
            {"visualProperty": "EDGE_TARGET_ARROW_SHAPE", "value": "NONE"},
        ],
        "mappings": [
            {
                "mappingType": "passthrough",
                "mappingColumn": "plot_label",
                "mappingColumnType": "String",
                "visualProperty": "NODE_LABEL",
            },
            {
                "mappingType": "discrete",
                "mappingColumn": "cluster_id",
                "mappingColumnType": "String",
                "visualProperty": "NODE_FILL_COLOR",
                "map": [{"key": key, "value": value} for key, value in cluster_colors.items()],
            },
            {
                "mappingType": "discrete",
                "mappingColumn": "pathway_source",
                "mappingColumnType": "String",
                "visualProperty": "NODE_SHAPE",
                "map": [
                    {"key": "Hallmark", "value": "DIAMOND"},
                    {"key": "Reactome", "value": "ELLIPSE"},
                ],
            },
            {
                "mappingType": "discrete",
                "mappingColumn": "representative_role",
                "mappingColumnType": "String",
                "visualProperty": "NODE_BORDER_WIDTH",
                "map": [
                    {"key": "Not selected", "value": "0.8"},
                    {"key": "Profile exemplar", "value": "4.5"},
                    {"key": "Central evidence", "value": "2.5"},
                ],
            },
            {
                "mappingType": "continuous",
                "mappingColumn": "negative_log10_minimum_FDR",
                "mappingColumnType": "Double",
                "visualProperty": "NODE_SIZE",
                "points": [
                    {"value": 0.0, "lesser": "18.0", "equal": "18.0", "greater": "18.0"},
                    {"value": 5.0, "lesser": "30.0", "equal": "30.0", "greater": "30.0"},
                    {"value": 20.0, "lesser": "48.0", "equal": "48.0", "greater": "48.0"},
                ],
            },
            {
                "mappingType": "continuous",
                "mappingColumn": "overlap_coefficient",
                "mappingColumnType": "Double",
                "visualProperty": "EDGE_WIDTH",
                "points": [
                    {"value": 0.50, "lesser": "0.5", "equal": "0.5", "greater": "0.5"},
                    {"value": 0.75, "lesser": "1.5", "equal": "1.5", "greater": "1.5"},
                    {"value": 1.00, "lesser": "3.0", "equal": "3.0", "greater": "3.0"},
                ],
            },
        ],
    }


def space_isolated_nodes(network_suid: int, view_suid: int, payload: dict[str, object]) -> None:
    degree = {node["data"]["id"]: 0 for node in payload["elements"]["nodes"]}
    for edge in payload["elements"]["edges"]:
        degree[edge["data"]["source"]] += 1
        degree[edge["data"]["target"]] += 1

    view = request(f"/networks/{network_suid}/views/{view_suid}")
    isolated = []
    connected = []
    for node in view["elements"]["nodes"]:
        pathway = node["data"]["name"]
        target = isolated if degree[pathway] == 0 else connected
        target.append(node)
    if not isolated or not connected:
        return

    x_values = [float(node["position"]["x"]) for node in connected]
    y_values = [float(node["position"]["y"]) for node in connected]
    x_target = max(x_values) + max(320.0, (max(x_values) - min(x_values)) * 0.28)
    y_min, y_max = min(y_values), max(y_values)
    isolated.sort(key=lambda node: (node["data"]["cluster_number"], node["data"]["name"]))
    usable_span = max(360.0, (y_max - y_min) * 0.70)
    step = usable_span / max(len(isolated) - 1, 1)
    y_start = (y_min + y_max - step * (len(isolated) - 1)) / 2

    for index, node in enumerate(isolated):
        properties = [
            {"visualProperty": "NODE_X_LOCATION", "value": x_target},
            {"visualProperty": "NODE_Y_LOCATION", "value": y_start + index * step},
        ]
        request(
            f"/networks/{network_suid}/views/{view_suid}/nodes/{node['data']['SUID']}?bypass=true",
            method="PUT",
            payload=properties,
        )


def main() -> None:
    version = request("/version")
    print(f"Cytoscape {version['cytoscapeVersion']} detected")
    FIGURES.mkdir(parents=True, exist_ok=True)
    command("session", "new")

    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    created = request("/networks?collection=COPD%20Schema%20pathway%20analysis", method="POST", payload=payload)
    network_suid = created["networkSUID"] if isinstance(created, dict) else created[0]
    views = request(f"/networks/{network_suid}/views")
    if not views:
        created_view = request(f"/networks/{network_suid}/views", method="POST")
        views = [created_view["networkViewSUID"]]
    view_suid = views[0]
    command("view", "set current", network=network_suid, view=view_suid)

    existing = request("/styles")
    if STYLE_NAME in existing:
        request(f"/styles/{quote(STYLE_NAME)}", method="DELETE")
    request("/styles", method="POST", payload=style_definition())
    time.sleep(1)
    request(f"/apply/styles/{quote(STYLE_NAME)}/{network_suid}")

    command(
        "layout",
        "force-directed",
        network=network_suid,
        edgeAttribute="overlap_coefficient",
        defaultSpringLength=75,
        defaultSpringCoefficient=0.000004,
        defaultNodeMass=3,
        numIterations=600,
        isDeterministic="true",
        singlePartition="false",
    )
    time.sleep(3)
    request(f"/apply/styles/{quote(STYLE_NAME)}/{network_suid}")
    time.sleep(1)
    space_isolated_nodes(network_suid, view_suid, payload)
    time.sleep(1)
    command("view", "fit content", view=view_suid)
    time.sleep(1)
    scale = request(
        f"/networks/{network_suid}/views/{view_suid}/network/NETWORK_SCALE_FACTOR"
    )["value"]
    request(
        f"/networks/{network_suid}/views/{view_suid}/network?bypass=true",
        method="PUT",
        payload=[{"visualProperty": "NETWORK_SCALE_FACTOR", "value": float(scale) * 0.82}],
    )
    time.sleep(1)

    stem = FIGURES / "FigS3_pathway_network_cytoscape"
    download(f"/networks/{network_suid}/views/{view_suid}.svg", stem.with_suffix(".svg"))
    download(f"/networks/{network_suid}/views/{view_suid}.pdf", stem.with_suffix(".pdf"))
    download(f"/networks/{network_suid}/views/{view_suid}.png?w=3600", stem.with_suffix(".png"))
    request(f"/session?{urlencode({'file': str(SESSION)})}", method="POST")
    print(f"Rendered network {network_suid}, view {view_suid} to {stem}.[svg|pdf|png]")
    print(f"Saved editable session: {SESSION}")


if __name__ == "__main__":
    main()

