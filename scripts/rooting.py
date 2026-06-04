#!/usr/bin/env python3
# Ported from OrthoFinder trees2ologs_of.py (Emms & Kelly 2019)
# S_IO/S_AD gene-tree rooting algorithm

import operator
import ete3


def _store_species_sets(tree, gene_map, tag="sp_"):
    tag_down = tag + "down"
    tag_up = tag + "up"
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature(tag_down, {gene_map(node.name)})
        elif node.is_root():
            continue
        else:
            node.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in node.get_children()]))
    for node in tree.traverse("preorder"):
        if node.is_root():
            node.add_feature(tag_up, set())
        else:
            parent = node.up
            others = [ch for ch in parent.get_children() if ch != node]
            sp_downs = set.union(*[other.__getattribute__(tag_down) for other in others])
            if parent.is_root():
                node.add_feature(tag_up, sp_downs)
            else:
                node.add_feature(tag_up, parent.__getattribute__(tag_up).union(sp_downs))
    tree.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in tree.get_children()]))


class _RootMap:
    def __init__(self, set_a, set_b, gene_to_species):
        self.set_a = set_a
        self.set_b = set_b
        self.gene_to_species = gene_to_species

    def gene_map(self, gene_name):
        sp = self.gene_to_species(gene_name)
        if sp in self.set_a:
            return True
        elif sp in self.set_b:
            return False
        else:
            raise ValueError(f"Species '{sp}' (from gene '{gene_name}') not in either clade set")


def _outgroup_ingroup_score(sp_up, sp_down, sett1, sett2, n_recip, n1, n2):
    f_dup = len(sp_up & sett1) * len(sp_up & sett2) * len(sp_down & sett1) * len(sp_down & sett2) * n_recip
    f_a = len(sp_up & sett1) * (n2 - len(sp_up & sett2)) * (n1 - len(sp_down & sett1)) * len(sp_down & sett2) * n_recip
    f_b = (n1 - len(sp_up & sett1)) * len(sp_up & sett2) * len(sp_down & sett1) * (n2 - len(sp_down & sett2)) * n_recip
    return max(f_dup, f_a, f_b)


def _get_roots(tree, species_tree_rooted, gene_to_species):
    """Return list of candidate root nodes (S_IO/S_AD algorithm)."""
    species_observed = {gene_to_species(g) for g in tree.get_leaf_names()}
    if len(species_observed) == 1:
        return [next(iter(tree.traverse()))]

    # Descend species tree to find deepest node with representatives in both children
    n = species_tree_rooted
    children = n.get_children()
    leaves = [set(ch.get_leaf_names()) for ch in children]
    have = [bool(l & species_observed) for l in leaves]
    while sum(have) < 2:
        n = children[have.index(True)]
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [bool(l & species_observed) for l in leaves]

    roots_list = []
    scores_list = []
    T, F, TF = {True}, {False}, {True, False}

    for i in range(len(leaves)):
        t1 = leaves[i]
        t2 = set.union(*[l for j, l in enumerate(leaves) if j != i])
        _store_species_sets(tree, gene_to_species)
        root_mapper = _RootMap(t1, t2, gene_to_species)
        sett1, sett2 = set(t1), set(t2)
        nt1, nt2 = float(len(t1)), float(len(t2))
        n_recip = 1.0 / (nt1 * nt1 * nt2 * nt2)
        _store_species_sets(tree, root_mapper.gene_map, "inout_")

        for m in tree.traverse("postorder"):
            if m.is_leaf():
                if len(m.inout_up) == 1 and m.inout_up != m.inout_down:
                    return [m]
            else:
                if len(m.inout_up) == 1 and len(m.inout_down) == 1 and m.inout_up != m.inout_down:
                    return [m]
                nodes = m.get_children() if m.is_root() else [m] + m.get_children()
                clades = [ch.inout_down for ch in nodes] if m.is_root() else ([m.inout_up] + [ch.inout_down for ch in m.get_children()])
                if len(nodes) == 3:
                    if all(len(c) == 1 for c in clades) and T in clades and F in clades:
                        if clades.count(T) == 1:
                            return [nodes[clades.index(T)]]
                        else:
                            return [nodes[clades.index(F)]]
                    elif T in clades and F in clades:
                        ab = [c == TF for c in clades]
                        idx = ab.index(True)
                        roots_list.append(nodes[idx])
                        scores_list.append(_outgroup_ingroup_score(nodes[idx].sp_up, nodes[idx].sp_down, sett1, sett2, n_recip, nt1, nt2))
                    elif clades.count(TF) >= 2:
                        roots_list.append(nodes[0])
                        scores_list.append(_outgroup_ingroup_score(nodes[0].sp_up, nodes[0].sp_down, sett1, sett2, n_recip, nt1, nt2))
                elif T in clades and F in clades:
                    roots_list.append(m)
                    scores_list.append(0)

    if not roots_list:
        return []
    return [sorted(zip(scores_list, roots_list), key=lambda x: x[0], reverse=True)[0][1]]


def root_gene_tree(tree, species_tree_rooted, gene_to_species):
    """Root an unrooted gene tree using the OrthoFinder S_IO/S_AD algorithm.

    Args:
        tree: ete3.Tree (unrooted, in-place modification)
        species_tree_rooted: ete3.Tree (rooted species tree, leaf names = species codes)
        gene_to_species: callable mapping leaf name -> species name/code

    Returns:
        True if rooting succeeded, False otherwise.
    """
    roots = _get_roots(tree, species_tree_rooted, gene_to_species)
    if not roots:
        return False
    root = max(roots, key=lambda r: r.get_closest_leaf()[1])
    if root is not tree:
        tree.set_outgroup(root)
    return True
