"""
Summarize 10,000 random networks by computing:
  1. Edge frequency  — for each undirected edge, count how many networks contained it,
                       then divide by 10,000.
  2. Average, min, max node degree — per node, across 10,000 networks.
  3. Average betweenness centrality (BC) — per node, across 10,000 networks.

Node degree and BC are read from each network's node_properties.txt.
Edge frequency is computed from each network's network.csv.

Usage:
    python summarize_networks.py

Outputs saved to:
    /nfs7/PHARM/Morgun_Lab/thanh/TkNA_analysis/random_nw/summary_output/
        edge_frequency.csv
        node_summary.csv
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

BASE_DIR   = Path('/nfs7/PHARM/Morgun_Lab/thanh/TkNA_analysis/random_nw/new_output_calc_bibc')
OUTPUT_DIR = Path('/nfs7/PHARM/Morgun_Lab/thanh/TkNA_analysis/random_nw/summary_output')
N_NETWORKS = 10000


# ── Parsers ───────────────────────────────────────────────────────────────────

def parse_node_properties(filepath):
    """
    Parse node_properties.txt which is TAB-separated.
    - Header line: 'name' [TAB] node1 [TAB] node2 ... (node names may contain spaces)
    - Data lines:  metric_name [TAB] value1 [TAB] value2 ...
    Returns: { metric_name: { node_name: float } }
    """
    with open(filepath, 'r') as f:
        lines = [line.rstrip('\n').rstrip('\t') for line in f
                 if line.strip() and not line.startswith('~')]
    if not lines:
        return {}

    # Split by tab — correctly handles node names with spaces
    header_tokens = lines[0].split('\t')
    node_names = header_tokens[1:]   # skip 'name' label

    data = {}
    for line in lines[1:]:
        tokens = line.split('\t')
        if not tokens:
            continue
        metric = tokens[0]
        values = tokens[1:]
        if len(values) == len(node_names):
            node_dict = {}
            for node, val in zip(node_names, values):
                try:
                    node_dict[node.strip()] = float(val.strip())
                except ValueError:
                    node_dict[node.strip()] = float('nan')
            data[metric] = node_dict
    return data


def parse_network_csv(filepath):
    """
    Parse network.csv and return a set of undirected edges as sorted tuples.
    (A, B) and (B, A) are treated as the same edge — undirected network.
    Auto-detects source/target columns; falls back to first two columns.
    """
    df   = pd.read_csv(filepath)
    cols = df.columns.tolist()

    src_col, tgt_col = None, None
    for c in cols:
        cl = c.lower()
        if cl in ('source', 'from', 'node1', 'src', 'cell1', 'celltype1'):
            src_col = c
        if cl in ('target', 'to', 'node2', 'tgt', 'cell2', 'celltype2'):
            tgt_col = c

    # Fallback: use first two columns
    if src_col is None or tgt_col is None:
        src_col, tgt_col = cols[0], cols[1]

    edges = set()
    for _, row in df.iterrows():
        u = str(row[src_col]).strip()
        v = str(row[tgt_col]).strip()
        if u != v:                              # skip self-loops
            edges.add(tuple(sorted([u, v])))   # sort so (A,B) == (B,A)
    return edges


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Accumulators (kept in memory across all 10,000 networks) ──────────────
    edge_counts = defaultdict(int)    # (node1, node2) -> count of networks containing this edge
    degree_sum  = defaultdict(float)  # node -> sum of degrees across networks
    degree_n    = defaultdict(int)    # node -> number of networks where node appeared
    degree_min  = {}                  # node -> minimum degree seen across networks
    degree_max  = {}                  # node -> maximum degree seen across networks
    bc_sum      = defaultdict(float)  # node -> sum of BC across networks
    bc_n        = defaultdict(int)    # node -> number of networks where node appeared

    n_processed = 0
    n_failed    = 0

    print(f"Processing {N_NETWORKS} networks from: {BASE_DIR}\n")

    for i in range(N_NETWORKS):
        folder    = BASE_DIR / str(i)
        net_file  = folder / 'network.csv'
        prop_file = folder / 'node_properties.txt'

        if not folder.exists():
            print(f"  [WARNING] Folder not found: {folder}")
            n_failed += 1
            continue

        # ── Edge frequency ─────────────────────────────────────────────────
        if net_file.exists():
            try:
                edges = parse_network_csv(net_file)
                for edge in edges:
                    edge_counts[edge] += 1    # each edge counted at most once per network
            except Exception as e:
                print(f"  [WARNING] Could not parse {net_file}: {e}")
        else:
            print(f"  [WARNING] Missing: {net_file}")

        # ── Node degree & BC ───────────────────────────────────────────────
        if prop_file.exists():
            try:
                props = parse_node_properties(prop_file)

                # Accumulate degree per node across networks
                for node, val in props.get('Node_degrees', {}).items():
                    if not np.isnan(val):
                        degree_sum[node] += val
                        degree_n[node]   += 1
                        # Track min and max degree per node
                        if node not in degree_min or val < degree_min[node]:
                            degree_min[node] = val
                        if node not in degree_max or val > degree_max[node]:
                            degree_max[node] = val

                # Accumulate BC per node across networks
                for node, val in props.get('Betweenness_centrality', {}).items():
                    if not np.isnan(val):
                        bc_sum[node] += val
                        bc_n[node]   += 1

            except Exception as e:
                print(f"  [WARNING] Could not parse {prop_file}: {e}")
        else:
            print(f"  [WARNING] Missing: {prop_file}")

        n_processed += 1
        if n_processed % 1000 == 0:
            print(f"  Processed {n_processed} / {N_NETWORKS} networks...")

    print(f"\nFinished. Networks processed: {n_processed}  |  Missing/failed: {n_failed}\n")

    # ── Output 1: Edge frequency ───────────────────────────────────────────────
    edge_rows = []
    for (u, v), count in sorted(edge_counts.items(), key=lambda x: -x[1]):
        edge_rows.append({
            'Node1':     u,
            'Node2':     v,
            'Count':     count,
            'Frequency': round(count / N_NETWORKS, 6)
        })

    edge_df  = pd.DataFrame(edge_rows)
    edge_out = OUTPUT_DIR / 'edge_frequency.csv'
    edge_df.to_csv(edge_out, index=False)
    print(f"[1] Edge frequency saved  →  {edge_out}  ({len(edge_df)} unique edges)")

    # ── Output 2: Average, min, max node degree & average BC per node ─────────
    all_nodes = sorted(set(list(degree_n.keys()) + list(bc_n.keys())))

    node_rows = []
    for node in all_nodes:
        avg_deg = (degree_sum[node] / degree_n[node]) if degree_n.get(node, 0) > 0 else np.nan
        min_deg = degree_min.get(node, np.nan)
        max_deg = degree_max.get(node, np.nan)
        avg_bc  = (bc_sum[node]    / bc_n[node])      if bc_n.get(node, 0)     > 0 else np.nan
        node_rows.append({
            'Node':                       node,
            'Avg_Node_Degree':            round(avg_deg, 6) if not np.isnan(avg_deg) else np.nan,
            'Min_Node_Degree':            int(min_deg) if not np.isnan(min_deg) else np.nan,
            'Max_Node_Degree':            int(max_deg) if not np.isnan(max_deg) else np.nan,
            'Avg_Betweenness_Centrality': round(avg_bc,  6) if not np.isnan(avg_bc)  else np.nan,
            'N_networks_degree':          degree_n.get(node, 0),
            'N_networks_BC':              bc_n.get(node, 0)
        })

    node_df  = pd.DataFrame(node_rows)
    node_out = OUTPUT_DIR / 'node_summary.csv'
    node_df.to_csv(node_out, index=False)
    print(f"[2] Node summary saved    →  {node_out}  ({len(node_df)} nodes)")

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
