import pandas as pd
import numpy as np
from collections import deque, defaultdict

DATA_PATH = "../pedigree/"

def build_parent_dict(df: pd.DataFrame) -> dict:
    ped = df.copy()
    ped["sampleID"] = ped["sampleID"].astype(str)
    ped["fatherID"] = ped["fatherID"].fillna(0)
    ped["motherID"] = ped["motherID"].fillna(0)

    def norm_parent(x):
        if pd.isna(x) or x == 0 or x == "0":
            return 0
        return str(x)

    ped["fatherID"] = ped["fatherID"].map(norm_parent)
    ped["motherID"] = ped["motherID"].map(norm_parent)

    parent_dict = {}
    for _, row in ped.iterrows():
        parents = []
        if row["fatherID"] != 0:
            parents.append(row["fatherID"])
        if row["motherID"] != 0:
            parents.append(row["motherID"])
        parent_dict[row["sampleID"]] = tuple(parents)
    return parent_dict


def add_dummy_founders(df: pd.DataFrame) -> pd.DataFrame:
    ped = df.copy()
    ped["sampleID"] = ped["sampleID"].astype(str)
    ped["fatherID"] = ped["fatherID"].fillna(0)
    ped["motherID"] = ped["motherID"].fillna(0)

    def norm_parent(x):
        if pd.isna(x) or x == 0 or x == "0":
            return 0
        return str(x)

    ped["fatherID"] = ped["fatherID"].map(norm_parent)
    ped["motherID"] = ped["motherID"].map(norm_parent)

    extra_rows = []
    dummy_counter = 1

    for idx, row in ped.iterrows():
        if row["fatherID"] == 0:
            dummy_id = f"__DUMMY_F_{dummy_counter}"
            dummy_counter += 1
            ped.at[idx, "fatherID"] = dummy_id
            extra_rows.append({"sampleID": dummy_id, "fatherID": 0, "motherID": 0})
        if row["motherID"] == 0:
            dummy_id = f"__DUMMY_M_{dummy_counter}"
            dummy_counter += 1
            ped.at[idx, "motherID"] = dummy_id
            extra_rows.append({"sampleID": dummy_id, "fatherID": 0, "motherID": 0})

    if extra_rows:
        ped = pd.concat([ped, pd.DataFrame(extra_rows)], ignore_index=True)

    ped = ped.drop_duplicates(subset=["sampleID"], keep="first").reset_index(drop=True)
    return ped


def topological_order(ped: pd.DataFrame) -> list:
    parent_map = {}
    for _, row in ped.iterrows():
        parents = [p for p in (row["fatherID"], row["motherID"]) if p != 0]
        parent_map[row["sampleID"]] = parents

    ordered = []
    remaining = set(parent_map.keys())

    while remaining:
        progress = False
        for sid in list(remaining):
            if all((p in ordered) or (p not in parent_map) for p in parent_map[sid]):
                ordered.append(sid)
                remaining.remove(sid)
                progress = True
        if not progress:
            raise ValueError("Pedigree contains a cycle or unresolved parent IDs.")
    return ordered


def kinship_matrix(df: pd.DataFrame):
    ped = add_dummy_founders(df)
    order = topological_order(ped)
    ped = ped.set_index("sampleID").loc[order].reset_index()

    ids = ped["sampleID"].tolist()
    id_to_idx = {sid: i for i, sid in enumerate(ids)}
    n = len(ids)

    phi = np.zeros((n, n), dtype=float)

    for i, row in ped.iterrows():
        father = row["fatherID"]
        mother = row["motherID"]

        is_founder = (father == 0 and mother == 0)
        if is_founder:
            phi[i, i] = 0.5
        else:
            f = id_to_idx[father]
            m = id_to_idx[mother]

            for j in range(i):
                phi[i, j] = 0.5 * (phi[f, j] + phi[m, j])
                phi[j, i] = phi[i, j]

            phi[i, i] = 0.5 * (1.0 + phi[f, m])

    return phi, id_to_idx


def ancestor_depths(person: str, parent_dict: dict) -> dict:
    """
    Return minimal generation distance from person to each ancestor.
    parent: distance 1, grandparent: distance 2, etc.
    """
    depths = {}
    q = deque([(person, 0)])
    visited = set()

    while q:
        node, d = q.popleft()
        if node in visited:
            continue
        visited.add(node)

        if d > 0:
            if node not in depths or d < depths[node]:
                depths[node] = d

        for p in parent_dict.get(node, []):
            q.append((p, d + 1))

    return depths


def exact_relationship(id1: str, id2: str, parent_dict: dict) -> str:
    parents1 = set(parent_dict.get(id1, ()))
    parents2 = set(parent_dict.get(id2, ()))

    # Direct ancestor / descendant
    anc1 = ancestor_depths(id1, parent_dict)
    anc2 = ancestor_depths(id2, parent_dict)

    if id2 in anc1:
        d = anc1[id2]
        if d == 1:
            return "parent-child"
        if d == 2:
            return "grandparent-grandchild"
        return f"{d-2}x great-grandparent / great-grandchild"

    if id1 in anc2:
        d = anc2[id1]
        if d == 1:
            return "parent-child"
        if d == 2:
            return "grandparent-grandchild"
        return f"{d-2}x great-grandparent / great-grandchild"

    # Siblings
    if len(parents1) == 2 and parents1 == parents2:
        return "full siblings"

    if len(parents1.intersection(parents2)) == 1:
        return "half siblings"

    # Avuncular
    # id1 is aunt/uncle of id2 if id1 is sibling of one of id2's parents
    for p in parents2:
        p_parents = set(parent_dict.get(p, ()))
        if len(parents1) == 2 and parents1 == p_parents:
            return "avuncular"
        if len(parents1.intersection(p_parents)) == 1:
            return "half-avuncular"

    for p in parents1:
        p_parents = set(parent_dict.get(p, ()))
        if len(parents2) == 2 and parents2 == p_parents:
            return "avuncular"
        if len(parents2.intersection(p_parents)) == 1:
            return "half-avuncular"

    # Cousins via nearest common ancestor depths
    common = set(anc1.keys()).intersection(anc2.keys())
    if common:
        best = min((max(anc1[a], anc2[a]), anc1[a], anc2[a], a) for a in common)
        _, d1, d2, _ = best

        # first cousins: both are grandchildren of common ancestor
        if d1 == 2 and d2 == 2:
            return "first cousins"
        if d1 == 2 and d2 == 3 or d1 == 3 and d2 == 2:
            return "first cousins once removed"
        if d1 == 3 and d2 == 3:
            return "second cousins"

        return f"common ancestor at depths ({d1}, {d2})"

    return "unrelated"


def pedigree_pairs(df: pd.DataFrame) -> pd.DataFrame:
    phi, id_to_idx = kinship_matrix(df)
    parent_dict = build_parent_dict(df)

    original_ids = df["sampleID"].astype(str).tolist()
    rows = []

    for i in range(len(original_ids)):
        for j in range(i + 1, len(original_ids)):
            id1 = original_ids[i]
            id2 = original_ids[j]
            kinship = phi[id_to_idx[id1], id_to_idx[id2]]
            pihat = 2.0 * kinship
            rel = exact_relationship(id1, id2, parent_dict)

            rows.append(
                {
                    "IID1": id1,
                    "IID2": id2,
                    "Relationship": rel,
                    "Kinship": kinship,
                    "Expected_PI_HAT": pihat,
                }
            )

    return pd.DataFrame(rows)


# Example
if __name__ == "__main__":
    df_genome = pd.read_csv("igsr-1000 genomes 30x on grch38-samples.tsv", sep='\t')
    pedigree = pd.read_csv("pedigree.txt", sep=" ")
    df_genome.columns = ["sampleID", "sex", "biosampleID", "pop", "pop_name", "subpop", "subpop_name", "popID", "collection"]
    df = df_genome.merge(pedigree, how = "left", on = "sampleID")
    pops = ["CHS", "YRI"]
    for pop in pops:
        df_pop = df[df["pop"] == pop].reset_index(drop=True).copy()
        out = pedigree_pairs(df_pop)
        out.to_csv(f"{DATA_PATH}pedigree_{pop}.csv", index=False)